'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import errno
import tarfile

import yaml
import pandas as pd
import pysam

import shutil

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


def generate_intervals(ref, size=100000000):
    fasta = pysam.FastaFile(ref)
    lengths = fasta.lengths
    names = fasta.references
 
    intervals = []
 
    for name, length in zip(names, lengths):
        if 'GL' in name:
            continue
        if 'MT' in name:
            continue

        for i in range((length/size)+1):
            intervals.append( name+ "_" + str(i*size) +"_"+ str((i+1)*size)) 
 

    return intervals

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
