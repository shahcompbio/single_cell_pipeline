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

    fastqs = pd.read_csv(fastqs_file, dtype=str)

    for column in ('cell_id', 'lane_id', 'fastq_1', 'fastq_2', 'source',):
        if column not in fastqs.columns:
            raise Exception(
                'input fastqs_file should contain {}'.format(column))

    if fastqs.duplicated(['cell_id', 'lane_id']).any():
        raise Exception('input fastqs_file with duplicate sample_id/lane_id pairs')

    fastq_1_filenames = dict()
    fastq_2_filenames = dict()
    for _, row in fastqs.iterrows():
        fastq_1_filenames[(row['cell_id'], row['lane_id'])] = row['fastq_1']
        fastq_2_filenames[(row['cell_id'], row['lane_id'])] = row['fastq_2']

    return fastq_1_filenames, fastq_2_filenames

def get_seqinfo(fastqs_file):

    fastqs = pd.read_csv(fastqs_file, dtype=str)

    for column in ('cell_id', 'source',):
        if column not in fastqs.columns:
            raise Exception(
                'input fastqs_file should contain {}'.format(column))

    seqinfo = {row["cell_id"]:row["source"] for _,row in fastqs.iterrows()}

    return seqinfo

def get_samples(fastqs_file):

    fastqs = pd.read_csv(fastqs_file, dtype=str)

    for column in ('cell_id',):
        if column not in fastqs.columns:
            raise Exception(
                'input fastqs_file should contain {}'.format(column))

    sample_ids = list(sorted(fastqs['cell_id'].unique()))
    return sample_ids


def get_bams(fastqs_file):

    fastqs = pd.read_csv(fastqs_file, dtype=str)

    for column in ('cell_id', 'bam',):
        if column not in fastqs.columns:
            raise Exception(
                'input bams_file should contain {}'.format(column))

    bam_filenames = dict()
    bai_filenames = dict()
    for _, row in fastqs.iterrows():
        bam_filenames[row['cell_id']] = row['bam']
        bai_filenames[row['cell_id']] = row['bam']+".bai"

    return bam_filenames, bai_filenames


def load_config(args):
    try:
        with open(args['config_file']) as infile:
            config = yaml.load(infile)

    except IOError:
        raise Exception(
            'Unable to open config file: {0}'.format(
                args['config_file']))
    return config



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
