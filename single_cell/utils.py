import yaml
import pandas as pd


def read_fastqs_file(args):

    fastqs = pd.read_csv(args['fastqs_file'], dtype=str)

    for column in ('sample_id', 'lane_id', 'fastq_1', 'fastq_2', 'source',):
        if column not in fastqs.columns:
            raise Exception(
                'input fastqs_file should contain {}'.format(column))

    sample_ids = list(sorted(fastqs['sample_id'].unique()))
    lanes = list(sorted(fastqs['lane_id'].unique()))

    fastq_1_filenames = dict()
    fastq_2_filenames = dict()
    for _, row in fastqs.iterrows():
        fastq_1_filenames[(row['sample_id'], row['lane_id'])] = row['fastq_1']
        fastq_2_filenames[(row['sample_id'], row['lane_id'])] = row['fastq_2']


    seqinfo = {row["sample_id"]:row["source"] for _,row in fastqs.iterrows()}

    return fastq_1_filenames, fastq_2_filenames, sample_ids, lanes, seqinfo


def load_config(args):
    try:
        with open(args['config_file']) as infile:
            config = yaml.load(infile)

    except IOError:
        raise Exception(
            'Unable to open config file: {0}'.format(
                args['config_file']))
    return config
