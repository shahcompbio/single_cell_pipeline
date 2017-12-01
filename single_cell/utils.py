import csv

import yaml
import pandas as pd
import pysam

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

def read_fastqs_file(fastqs_file):

    fastqs = pd.read_csv(fastqs_file, dtype=str)

    for column in ('sample_id', 'lane_id', 'fastq_1', 'fastq_2', 'source',):
        if column not in fastqs.columns:
            raise Exception(
                'input fastqs_file should contain {}'.format(column))

    sample_ids = list(sorted(fastqs['sample_id'].unique()))

    if fastqs.duplicated(['sample_id', 'lane_id']).any():
        raise Exception('input fastqs_file with duplicate sample_id/lane_id pairs')

    fastq_1_filenames = dict()
    fastq_2_filenames = dict()
    for _, row in fastqs.iterrows():
        fastq_1_filenames[(row['sample_id'], row['lane_id'])] = row['fastq_1']
        fastq_2_filenames[(row['sample_id'], row['lane_id'])] = row['fastq_2']

    seqinfo = {row["sample_id"]:row["source"] for _,row in fastqs.iterrows()}

    return fastq_1_filenames, fastq_2_filenames, sample_ids, seqinfo


def load_config(args):
    try:
        with open(args['config_file']) as infile:
            config = yaml.load(infile)

    except IOError:
        raise Exception(
            'Unable to open config file: {0}'.format(
                args['config_file']))
    return config


def concatenate_csv(in_filenames, out_filename):
    """merge csv files, uses csv module to handle inconsistencies in column
    indexes, pandas uses a lot of memory
    :param in_filenames: input file dict
    :param out_filename: output file
    """
    writer = None
    for _,infile in in_filenames.iteritems():

        with open(infile) as inp:
            reader= csv.DictReader(inp)

            for row in reader:
                if not writer:
                    writer = csv.DictWriter(open(out_filename, "w"),
                                            fieldnames=reader._fieldnames)
                    writer.writeheader()

                writer.writerow(row)
