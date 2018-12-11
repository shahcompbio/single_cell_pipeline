'''
Created on Jun 6, 2018

@author: dgrewal
'''
import single_cell
import os
import yaml
import collections

import pipeline_config


class folded_unicode(unicode):
    pass


class literal_unicode(unicode):
    pass



def folded_unicode_representer(dumper, data):
    return dumper.represent_scalar(u'tag:yaml.org,2002:str', data, style='>')


def literal_unicode_representer(dumper, data):
    return dumper.represent_scalar(u'tag:yaml.org,2002:str', data, style='|')


yaml.add_representer(folded_unicode, folded_unicode_representer)
yaml.add_representer(literal_unicode, literal_unicode_representer)


def override_config(config, override):
    def update(d, u):
        for k, v in u.iteritems():
            if isinstance(v, collections.Mapping):
                d[k] = update(d.get(k, {}), v)
            else:
                d[k] = v
        return d

    if not override:
        return config

    cfg = update(config, override)

    return cfg


def get_batch_params(override=None):

    pools = {
        "standard": {
            "tasks_per_node": 4,
            "numcores": 1,
            "mem": 8,
            "dedicated": False
        },
        "highmem": {
            "tasks_per_node": 2,
            "numcores": 1,
            "mem": 16,
            "dedicated": False
        },
        "multicore": {
            "tasks_per_node": 1,
            "numcores": 8,
            "mem": 8,
            "dedicated": False
        }
    }

    data = {
        "storage_container_name": "tasks-container",
        "no_delete_pool": True,
        "no_delete_job": False,
        "pools": pools,
        "reference": "grch37"
    }

    data = override_config(data, override)

    return data


def generate_autoscale_formula(tasks_per_node, dedicated):

    if dedicated:
        node_type = "TargetDedicatedNodes"
    else:
        node_type = "TargetLowPriorityNodes"

    formula = (
        "tasksPerNode = {0};\n"
        "numAddMax = 20;\n"
        "numDelMax = 20;\n"
        "startingNumberOfVMs = 0;\n"
        "minNumberofVMs = 0;\n"
        "maxNumberofVMs = 1000;\n"
        "pendingTaskSamplePercent = $PendingTasks.GetSamplePercent(180 * TimeInterval_Second);\n"
        "pendingTaskSamples = pendingTaskSamplePercent < 70 ? startingNumberOfVMs : avg($PendingTasks.GetSample(180 * TimeInterval_Second));\n"
        "cores = ${1} * tasksPerNode;\n"
        "$extraVMs = (pendingTaskSamples - cores) / tasksPerNode;\n"
        "$extraVMs = $extraVMs + (tasksPerNode-1)/tasksPerNode;\n"
        "$extraVMs = min(numAddMax, $extraVMs);\n"
        "$extraVMs = max(-numDelMax, $extraVMs);\n"
        "targetVMs = (${1} + $extraVMs);\n"
        "${1} = max(minNumberofVMs,min(maxNumberofVMs, targetVMs));\n"
    )

    formula = formula.format(tasks_per_node, node_type)

    formula = literal_unicode(formula)

    return formula


def create_vm_commands():

    # GETSIZE AND MOUNT:
    # the disks are not guaranteed to appear under the same device name when starting a new
    # vm from a managed image. So we cannot be sure where the devices should be mounted,
    # this code uses the disk size to guess where to mount
    # If disk in under 40ish GB then mount it to refdata (our default reference disk is 32 GB)
    # If disk is over 900ish GB then mount it to datadrive (our default scratch space is 1TB)
    # DOCKER GROUP:
    # We do not have sudo available in tasks, docker requires sudo. Docker also comes with a
    # 'docker' user group. This user group has escalated permissions, any user in this group
    # can run docker without sudo (but still as an escalated user). Add the current user (default _azbatch)
    # to the docker user group. This cannot be done during the image creation, image capture removes all user
    # information.
    commands = (
        "sudo gpasswd -a $USER docker\n"
        "sudo mkdir /datadrive\n"
        "sudo chmod -R 777 /datadrive\n"
    )

    commands = literal_unicode(commands)

    return commands


def get_vm_size_azure(numcores):
    if numcores <= 2:
        return "STANDARD_E2_V3"
    elif numcores <= 4:
        return "STANDARD_E4_V3"
    elif numcores <= 8:
        return "STANDARD_E8_V3"
    else:
        # max 16 cores
        return "STANDARD_E16_V3"


def get_vm_image_id():
    subscription = os.environ.get("SUBSCRIPTION_ID", "id-missing")
    resource_group = os.environ.get("RESOURCE_GROUP", "id-missing")
    imageid = ['subscriptions', subscription, 'resourceGroups',
               resource_group, 'providers', 'Microsoft.Compute',
               'images', 'dockerproduction-smalldisk']
    imageid = '/'.join(imageid)
    imageid = '/' + imageid
    return imageid

def get_pool_def(
        tasks_per_node, reference, pool_type, numcores, memory, dedicated):

    autoscale_formula = generate_autoscale_formula(tasks_per_node, dedicated)

    vm_commands = create_vm_commands()

    poolname = "singlecell{}{}".format(reference, pool_type)

    pooldata = {
        "pool_vm_size": get_vm_size_azure(tasks_per_node),
        "cpus": numcores,
        'mem': memory,
        'dedicated': dedicated,
        'node_resource_id': get_vm_image_id(),
        'node_os_publisher': 'Canonical',
        'node_os_offer': 'UbuntuServer',
        'node_os_sku': 'batch.node.ubuntu 16.04',
        'data_disk_sizes': None,
        'max_tasks_per_node': tasks_per_node,
        'auto_scale_formula': autoscale_formula,
        'create_vm_commands': vm_commands,
        'start_resources': None
    }

    return {poolname: pooldata}


def get_compute_start_commands():
    # Azure batch starts a process on each node during startup
    # this process runs the entire time VM is up (speculation?). running sudo gasswd
    # adds our user to the docker group in node startup. but all user/group
    # permission updates require restarting the bash process (normally by
    # running newgrp or by logout-login). Since we cannot do that with batch
    # and since the start task and compute task run under same process, we need
    # to explicitly specify the user group when running docker commands. sg in
    # linux can be used to set group when executing commands

    commands = (
        'clean_up () {',
        '  echo "clean_up task executed"',
        '  find $AZ_BATCH_TASK_WORKING_DIR/ -xtype l -delete',
        '  exit 0',
        '}',
        'trap clean_up EXIT',
        'mkdir -p /datadrive/$AZ_BATCH_TASK_WORKING_DIR/\n'
        'cd /datadrive/$AZ_BATCH_TASK_WORKING_DIR/\n'
    )

    commands = '\n'.join(commands)

    commands = literal_unicode(commands)

    return {"compute_start_commands": commands}


def get_compute_finish_commands():
    commands = (
        'sg docker -c "docker run -v /datadrive:/datadrive -w /datadrive/$AZ_BATCH_TASK_WORKING_DIR ubuntu find . -xtype l -delete"\n'
        'sg docker -c "docker run -v /datadrive:/datadrive -w /datadrive/$AZ_BATCH_TASK_WORKING_DIR ubuntu find . -xtype f -delete"\n'
    )

    commands = literal_unicode(commands)

    return {"compute_finish_commands": commands}


def get_all_pools(pool_config, reference):

    pooldefs = {}

    for pooltype, poolinfo in pool_config.iteritems():
        tasks_per_node = poolinfo["tasks_per_node"]
        numcores = poolinfo["numcores"]
        memory = poolinfo["mem"]
        dedicated = poolinfo["dedicated"]

        pool_def = get_pool_def(
            tasks_per_node, reference, pooltype,
            numcores, memory, dedicated
        )

        pooldefs.update(pool_def)

    return {"pools": pooldefs}


def write_config(params, filepath):
    with open(filepath, 'w') as outputfile:
        yaml.dump(params, outputfile, default_flow_style=False)


def get_batch_config(defaults, override=None):
    config = {}

    config.update(
        get_all_pools(
            defaults['pools'],
            defaults['reference']))

    config.update(
        {"storage_container_name": defaults["storage_container_name"]}
    )

    config.update(get_compute_start_commands())
    config.update(get_compute_finish_commands())

    config.update({"no_delete_pool": defaults["no_delete_pool"]})
    config.update({"no_delete_job": defaults["no_delete_job"]})

    config = override_config(config, override)

    return config
