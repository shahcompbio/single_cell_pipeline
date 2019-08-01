import os
import single_cell
import subprocess
import collections
import numbers
import yaml


def walk_yaml(yamldata):
    mounts = []

    for k, v in yamldata.items():
        if isinstance(v, collections.Mapping):
            mounts += walk_yaml(v)

        if isinstance(v, basestring) or isinstance(v, numbers.Number):
            mounts.append(v)

        if isinstance(v, collections.Sequence) and not isinstance(v, basestring):
            for val in v:
                mounts += walk_yaml(val)

    return mounts


def load_yaml(path):
    try:
        with open(path) as infile:
            data = yaml.safe_load(infile)

    except IOError:
        raise Exception(
            'Unable to open file: {0}'.format(path))
    return data


def get_docker():
    for path in os.environ["PATH"].split(os.pathsep):
        if os.path.exists(os.path.join(path, 'docker')) and os.path.isfile(os.path.join(path, 'docker')):
            return os.path.join(path, 'docker')
    return None


def is_path(path_str):
    if not isinstance(path_str, basestring):
        return

    if '/' in path_str:
        return True
    # if its a directory or file in current dir
    if os.path.exists(path_str):
        return True


def is_yaml(path_str):
    extension = os.path.splitext(path_str)[-1]

    if extension in ['.yaml', '.yml']:
        return True


def get_root_dir(path_str):
    path_str = os.path.abspath(path_str)
    path_split = path_str.split('/')

    return '/' + path_split[1]


def get_mounts_from_yaml(yamlfile):
    data = load_yaml(yamlfile)

    mounts = []

    leaf_nodes = walk_yaml(data)

    for leafnode in leaf_nodes:
        if is_path(leafnode):
            mounts.append(get_root_dir(leafnode))

    return mounts


def get_path_mounts(args):
    mounts = set()

    for arg, argvalue in args.items():
        if is_path(argvalue):
            mounts.add(get_root_dir(argvalue))
            if is_yaml(argvalue):
                yamlmounts = get_mounts_from_yaml(argvalue)
                for mnt in yamlmounts:
                    mounts.add(mnt)

    return list(mounts)


def get_mounts(args):
    mounts = []
    docker_path = get_docker()

    docker_socket = '/var/run/docker.sock'

    mounts.append(docker_path)
    mounts.append(docker_socket)

    mounts.extend(get_path_mounts(args))

    return mounts


def get_env_vars():
    env_vars = [
        'AZURE_BATCH_URL',
        'CLIENT_ID',
        'TENANT_ID',
        'SECRET_KEY',
        'SUBSCRIPTION_ID',
        'RESOURCE_GROUP',
        'AZURE_KEYVAULT_ACCOUNT',
    ]

    return env_vars


def get_version():
    version = single_cell.__version__
    version = version.split('+')
    return version[0]


def get_docker_command(args):
    workingdir = os.getcwd()

    mounts = get_mounts(args)

    env_vars = get_env_vars()

    version = get_version()

    container = "shahlab.azurecr.io/scp/single_cell_pipeline:{}".format(version)

    docker_path = get_docker()

    cmd = [
        docker_path, 'run', '-w', workingdir, '--rm'
    ]

    for mount in mounts:
        cmd.extend(['-v', '{0}:{0}'.format(mount)])

    for env_var in env_vars:
        cmd.extend(['-e', env_var])

    cmd.append(container)

    return cmd


def run_with_docker(args, sc_cmd):
    cmd = get_docker_command(args)

    cmd.extend(sc_cmd)

    docker_part = cmd.index('--run_with_docker')
    del cmd[docker_part]

    print ("running command: {}".format(' '.join(cmd)))

    subprocess.check_call(cmd)
