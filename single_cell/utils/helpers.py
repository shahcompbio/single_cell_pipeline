'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import errno
import tarfile

import yaml
import pysam

import shutil

from collections import OrderedDict
import warnings
import subprocess
from subprocess import Popen, PIPE

import multiprocessing


from multiprocessing.pool import ThreadPool

from single_cell.config import generate_pipeline_config
from single_cell.config import generate_batch_config

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


def generate_configs_in_temp(args):

    config_yaml = "config.yaml"
    batch_yaml = "batch.yaml"

    pipelineargs = args[args.keys()[0]]

    tmpdir = pipelineargs.get("tmpdir", None)
    #use pypeliner tmpdir to store yaml
    if tmpdir:
        config_yaml = os.path.join(tmpdir, config_yaml)
        batch_yaml = os.path.join(tmpdir, batch_yaml)
    else:
        warnings.warn("no tmpdir specified, generating configs in working dir")
        config_yaml = os.path.join(os.getcwd(), config_yaml)
        batch_yaml = os.path.join(os.getcwd(), batch_yaml)

    config_yaml = get_incrementing_filename(config_yaml)
    batch_yaml = get_incrementing_filename(batch_yaml)

    params_override = pipelineargs["config_override"]

    generate_pipeline_config.main(output=config_yaml, input_params = params_override)

    generate_batch_config.main(output=batch_yaml, input_params = params_override)

    for _, mode_args in args.iteritems():
        mode_args["config_file"] = config_yaml
        mode_args["submit_config"] = batch_yaml

    return args

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

    stdout=PIPE
    if output:
        stdout=open(output, "w")

    p = Popen(cmd, stdout=stdout, stderr=PIPE)

    cmdout, cmderr = p.communicate()
    retc  = p.returncode

    if retc:
        raise Exception("command failed. stderr:{}, stdout:{}".format(cmdout, cmderr))

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


def get_regions(chromosome_lengths, split_size):
    if split_size is None:
        return dict(enumerate(chromosome_lengths.keys()))

    regions = []

    for chrom, length in chromosome_lengths.iteritems():
        lside_interval = range(1, length + 1, split_size)
        rside_interval = range(split_size, length + split_size, split_size)

        for beg, end in zip(lside_interval, rside_interval):
            end = min(end, length)

            regions.append('{}-{}-{}'.format(chrom, beg, end))

    return regions


def load_bam_chromosome_lengths(file_name, chromosomes=None):

    chromosome_lengths = OrderedDict()

    bam = pysam.Fastafile(file_name)

    for chrom, length in zip(bam.references, bam.lengths):
        if chromosomes and chrom not in chromosomes:
            continue

        chromosome_lengths[str(chrom)] = int(length)

    return chromosome_lengths


def get_bam_regions(bam_file, split_size, chromosomes):
    chromosome_lengths = load_bam_chromosome_lengths(bam_file, chromosomes=chromosomes)
    return get_regions(chromosome_lengths, split_size)



def get_fastqs(fastqs_file):

    data = load_yaml(fastqs_file)

    for cell in data.keys():
        assert "fastqs" in data[cell], "couldnt extract fastq file paths from yaml input for cell: {}".format(cell)

    fastq_1_filenames = dict()
    fastq_2_filenames = dict()
    for cell in data.keys():
        fastqs = data[cell]["fastqs"]
        
        for lane,paths in fastqs.iteritems():
            fastq_1_filenames[(cell,lane)] = paths["fastq_1"]
            fastq_2_filenames[(cell,lane)] = paths["fastq_2"]
    
    return fastq_1_filenames, fastq_2_filenames

def get_seqinfo(fastqs_file):

    data = load_yaml(fastqs_file)

    for cell in data.keys():
        assert "fastqs" in data[cell], "couldnt extract fastq file paths from yaml input for cell: {}".format(cell)

    seqinfo = dict()
    for cell in data.keys():
        fastqs = data[cell]["fastqs"]


        for lane,paths in fastqs.iteritems():

            if "source" not in paths:
                raise Exception("source key missing in cell: {}".format(cell))
            seqinfo[(cell,lane)] = paths["source"]

    return seqinfo


def get_sample_info(fastqs_file):
    """
    load yaml and remove some extra info to reduce size
    """

    data = load_yaml(fastqs_file)

    cells = data.keys()
    
    for cell in cells:
        del data[cell]["fastqs"]
        del data[cell]["bam"]
    return data


def get_samples(fastqs_file):

    data = load_yaml(fastqs_file)
    
    return data.keys()

def get_bams(fastqs_file):

    data = load_yaml(fastqs_file)

    for cell in data.keys():
        assert "bam" in data[cell], "couldnt extract bam file paths from yaml input for cell: {}".format(cell)

    bam_filenames = {cell:data[cell]["bam"] for cell in data.keys()}
    bai_filenames = {cell:data[cell]["bam"]+".bai" for cell in data.keys()}

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
