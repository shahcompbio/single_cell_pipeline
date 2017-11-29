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

def get_seq_info(fastqs_file, sample_id, lane):
    fastqs = pd.read_csv(fastqs_file, dtype=str)

    for column in ('sample_id', 'lane_id', 'fastq_1', 'fastq_2', 'source',):
        if column not in fastqs.columns:
            raise Exception(
                'input fastqs_file should contain {}'.format(column))

    return fastqs.loc[(fastqs['sample_id'] == sample_id) & (fastqs['lane_id'] == lane)].iloc[0]["source"]


def read_fastqs_file(args):

    fastqs = pd.read_csv(args['fastqs_file'], dtype=str)

    for column in ('sample_id', 'lane_id', 'fastq_1', 'fastq_2', 'source',):
        if column not in fastqs.columns:
            raise Exception(
                'input fastqs_file should contain {}'.format(column))

    sample_ids = list(sorted(fastqs['sample_id'].unique()))

    fastq_1_filenames = dict()
    fastq_2_filenames = dict()
    for _, row in fastqs.iterrows():
        fastq_1_filenames[(row['sample_id'], row['lane_id'])] = row['fastq_1']
        fastq_2_filenames[(row['sample_id'], row['lane_id'])] = row['fastq_2']



    return fastq_1_filenames, fastq_2_filenames, sample_ids


def load_config(args):
    try:
        with open(args['config_file']) as infile:
            config = yaml.load(infile)

    except IOError:
        raise Exception(
            'Unable to open config file: {0}'.format(
                args['config_file']))
    return config
