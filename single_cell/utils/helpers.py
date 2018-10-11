'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import errno
import tarfile
import time
import yaml

import shutil

from subprocess import Popen, PIPE

import multiprocessing

from multiprocessing.pool import ThreadPool


def get_coltype_reference():
    coltypes = {
        'estimated_library_size': 'int', 'total_mapped_reads': 'int',
        'total_reads_hmmcopy': 'float', 'cor_map': 'float',
        'cv_neutral_state': 'float', 'MBRSM_dispersion': 'float',
        'map': 'float', 'mad_hmmcopy': 'float', 'copy': 'float',
        'modal_curve': 'float', 'sample_well': 'str', 'true_multiplier': 'float',
        'reads': 'int', 'jira_id': 'str', 'gc': 'float', 'integer_copy_number': 'int',
        'breakpoints': 'int', 'total_duplicate_reads': 'int',
        'quality': 'float', 'cv_hmmcopy': 'float', 'empty_bins_hmmcopy': 'float',
        'paired_mapped_reads': 'int', 'total_reads': 'int', 'end': 'int', 'width': 'int',
        'MBRSI_dispersion_non_integerness': 'float', 'total_properly_paired': 'int',
        'chr': 'str', 'sample_type': 'str', 'mean_insert_size': 'float', 'start': 'int',
        'state': 'int', 'valid': 'bool', 'coverage_breadth': 'float', 'empty_bins_hmmcopy_chrY': 'int',
        'too_even': 'bool', 'unpaired_duplicate_reads': 'int', 'unpaired_mapped_reads': 'int',
        'unmapped_reads': 'int', 'mad_chr19': 'float', 'cell_id': 'str', 'cell_call': 'str',
        'coverage_depth': 'float', 'median_insert_size': 'float', 'modal_quantile': 'float',
        'sample_plate': 'str', 'mean_state_mads': 'float', 'ideal': 'bool',
        'experimental_condition': 'str', 'mean_copy': 'float', 'mean_hmmcopy_reads_per_bin': 'float',
        'multiplier': 'int', 'percent_duplicate_reads': 'float', 'i7_barcode': 'str',
        'total_halfiness': 'float', 'std_hmmcopy_reads_per_bin': 'float',
        'standard_deviation_insert_size': 'float', 'mean_state_vars': 'float',
        'all_heatmap_order': 'int', 'scaled_halfiness': 'float', 'cor_gc': 'float',
        'median': 'float', 'state_mode': 'int', 'paired_duplicate_reads': 'int',
        'median_hmmcopy_reads_per_bin': 'float', 'mad_neutral_state': 'float',
        'autocorrelation_hmmcopy': 'float', 'mad_autosomes': 'float', 'i5_barcode': 'str',
        'loglikehood': 'float', 'MSRSI_non_integerness': 'float'}

    ignore_cols = set(range(9))

    return coltypes, ignore_cols


def resolve_template(regions, template, format_key):
    outputs = {v: template.format(**{format_key:v}) for v in regions}
    return outputs

def get_fastq_files(input_yaml):
    data = load_yaml(input_yaml)

    items = {}
    for cell_id, cell_data in data.iteritems():
        items[cell_id] = {}
        for lane, laneinfo in cell_data["fastqs"].iteritems():
            items[cell_id][lane] = {}
            items[cell_id][lane]['fastq_1'] = format_file_yaml(laneinfo['fastq_1'])
            items[cell_id][lane]['fastq_2'] = format_file_yaml(laneinfo['fastq_2'])
    return items

def format_file_yaml(filepath):
    ext = os.path.splitext(filepath)

    if ext[1] == ".gz":
        ext = os.path.splitext(ext[0])

    mapping = {'.bam':'bam', '.pdf': 'PDF',
               '.fastq': 'fastq', '.h5': 'H5',
               '.tar': 'tar'}

    return {'filename': filepath, 'type': mapping[ext[1]]}

def get_container_ctx(container_config, image_name, docker_only=False):
    if docker_only and not container_config['container_type'] == 'docker':
        return {}

    credentials = container_config['images'][image_name]
    docker_context = {
              'image': credentials['image'],
              'container_type': container_config['container_type'],
              'mounts': container_config['mounts'],
              'username': credentials['username'],
              'password': credentials['password'],
              'server': credentials['server'],
          }
    return docker_context


def get_mount_dirs_docker(*args):
    mounts = set()

    for arg in args:
        if os.path.exists(os.path.dirname(arg)):
            if not arg.startswith('/'):
                arg = os.path.abspath(arg)
            arg = arg.split('/')

            mounts.add('/' + arg[1])
    return sorted(mounts)


def write_to_yaml(outfile, data):
    with open(outfile, 'w') as output:
        yaml.safe_dump(data, output, default_flow_style=False)


def eval_expr(val, operation, threshold):

    if operation == "gt":
        if val > threshold:
            return True
    elif operation == 'ge':
        if val >= threshold:
            return True
    elif operation == 'lt':
        if val < threshold:
            return True
    elif operation == 'le':
        if val <= threshold:
            return True
    elif operation == 'eq':
        if val == threshold:
            return True
    elif operation == 'ne':
        if not val == threshold:
            return True
    else:
        raise Exception("unknown operator type: {}".format(operation))

    return False


def get_incrementing_filename(path):
    """
    avoid overwriting files. if path exists then return path
    otherwise generate a path that doesnt exist.
    """

    if not os.path.exists(path):
        return path

    i = 0
    while os.path.exists("{}.{}".format(path, i)):
        i += 1

    return "{}.{}".format(path, i)


def run_in_parallel(worker, args, ncores=None):

    def args_unpack(worker, args):
        return worker(*args)

    count = multiprocessing.cpu_count()

    if ncores:
        count = min(ncores, count)

    pool = ThreadPool(processes=count)

    tasks = []

    for arg in args:

        task = pool.apply_async(args_unpack,
                                args=(worker, arg),
                                )
        tasks.append(task)

    pool.close()
    pool.join()

    [task.get() for task in tasks]

    pool.terminate()
    del pool


def run_cmd(cmd, output=None):

    stdout = PIPE
    if output:
        stdout = open(output, "w")

    p = Popen(cmd, stdout=stdout, stderr=PIPE)

    cmdout, cmderr = p.communicate()
    retc = p.returncode

    if retc:
        raise Exception(
            "command failed. stderr:{}, stdout:{}".format(
                cmdout,
                cmderr))

    if output:
        stdout.close()


def load_yaml(path):
    try:
        with open(path) as infile:
            data = yaml.load(infile)

    except IOError:
        raise Exception(
            'Unable to open file: {0}'.format(path))
    return data


def symlink(actual_file, symlink):
    if not os.path.exists(symlink):
        os.symlink(actual_file, symlink)


def copy_file(infile, output):
    shutil.copy(infile, output)


def get_fastqs(fastqs_file):

    data = load_yaml(fastqs_file)

    for cell in data.keys():
        assert "fastqs" in data[
            cell], "couldnt extract fastq file paths from yaml input for cell: {}".format(cell)

    fastq_1_filenames = dict()
    fastq_2_filenames = dict()
    for cell in data.keys():
        fastqs = data[cell]["fastqs"]

        for lane, paths in fastqs.iteritems():
            fastq_1_filenames[(cell, lane)] = paths["fastq_1"]
            fastq_2_filenames[(cell, lane)] = paths["fastq_2"]

    return fastq_1_filenames, fastq_2_filenames


def get_instrument_info(fastqs_file):

    data = load_yaml(fastqs_file)

    for cell in data.keys():
        assert "fastqs" in data[
            cell], "couldnt extract fastq file paths from yaml input for cell: {}".format(cell)

    seqinfo = dict()
    for cell in data.keys():
        fastqs = data[cell]["fastqs"]

        for lane, paths in fastqs.iteritems():

            if "sequencing_instrument" not in paths:
                raise Exception(
                    "instrument key missing in cell: {}".format(cell))
            seqinfo[(cell, lane)] = paths["sequencing_instrument"]

    return seqinfo


def get_center_info(fastqs_file):

    data = load_yaml(fastqs_file)

    for cell in data.keys():
        assert "fastqs" in data[
            cell], "couldnt extract fastq file paths from yaml input for cell: {}".format(cell)

    seqinfo = dict()
    for cell in data.keys():
        fastqs = data[cell]["fastqs"]

        for lane, paths in fastqs.iteritems():

            if "sequencing_center" not in paths:
                raise Exception(
                    "sequencing_center key missing in cell: {}".format(cell))
            seqinfo[(cell, lane)] = paths["sequencing_center"]

    return seqinfo


def get_sample_info(fastqs_file):
    """
    load yaml and remove some extra info to reduce size
    """

    data = load_yaml(fastqs_file)

    cells = data.keys()

    for cell in cells:
        data[cell]["cell_call"] = data[cell]["pick_met"]
        data[cell]["experimental_condition"] = data[cell]["condition"]
        del data[cell]["fastqs"]
        del data[cell]["bam"]
        del data[cell]["pick_met"]
        del data[cell]["condition"]
    
    return data


def get_samples(fastqs_file):

    data = load_yaml(fastqs_file)

    return data.keys()


def get_bams(fastqs_file):

    data = load_yaml(fastqs_file)

    for cell in data.keys():
        assert "bam" in data[
            cell], "couldnt extract bam file paths from yaml input for cell: {}".format(cell)

    bam_filenames = {cell: data[cell]["bam"] for cell in data.keys()}
    bai_filenames = {cell: data[cell]["bam"] + ".bai" for cell in data.keys()}

    return bam_filenames, bai_filenames


def load_config(args):

    return load_yaml(args["config_file"])


def makedirs(directory, isfile=False):

    if isfile:
        directory = os.path.dirname(directory)

    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))
