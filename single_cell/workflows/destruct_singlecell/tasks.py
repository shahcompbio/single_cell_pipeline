import gzip
import os
import random
from itertools import islice

import pandas as pd
import pypeliner


def destruct_bamdisc(
        destruct_config, normal_bam_file, stats,
        reads_1, reads_2, sample_1, sample_2, tempdir, cell_id,
        tag=False
):
    """
    :param destruct_config: dict object with destruct params
    :type destruct_config: dict
    :param normal_bam_file: normal bam file path
    :type normal_bam_file: str
    :param stats: path to output stats file
    :type stats: str
    :param reads_1: path to output R1 fastq
    :type reads_1: str
    :param reads_2: path to output R2 fastq
    :type reads_2: str
    :param sample_1: path to output R1 sample fastq
    :type sample_1: str
    :param sample_2: path to output R2 sample fastq
    :type sample_2: str
    :param tempdir: path to dir for storing temp file
    :type tempdir: str
    :param cell_id: cell_id to tag reads with
    :type cell_id: str
    :param tag: Tag reads in fastq with cell_id if set
    :type tag: bool
    """

    pypeliner.helpers.makedirs(tempdir)

    if tag:
        disc_reads_1 = os.path.join(tempdir, "pre_tagged_r1.fastq.gz")
        disc_reads_2 = os.path.join(tempdir, "pre_tagged_r2.fastq.gz")
    else:
        disc_reads_1 = reads_1
        disc_reads_2 = reads_2

    cmd = ['destruct_bamdiscordantfastq',
           '-r',
           '-c', destruct_config['bam_max_soft_clipped'],
           '-f', destruct_config['bam_max_fragment_length'],
           '-b', normal_bam_file,
           '-s', stats,
           '--fastq1', disc_reads_1,
           '--fastq2', disc_reads_2,
           '-t', tempdir,
           '-n', destruct_config['num_read_samples'],
           '--sample1', sample_1,
           '--sample2', sample_2,
           ]
    pypeliner.commandline.execute(*cmd)

    if tag:
        tag_reads(disc_reads_1, reads_1, cell_id)
        tag_reads(disc_reads_2, reads_2, cell_id)


def tag_reads(input_fastq, output_fastq, tag):
    """
    tag the reads in the fastq file with the provided tag
    :param input_fastq: path to fastq file
    :type input_fastq: str
    :param output_fastq: path to output tagged fastq
    :type output_fastq: str
    :param tag: value to tag reads with
    :type tag: str
    """
    opener = (open, gzip.open)[input_fastq.endswith('.gz')]
    with opener(input_fastq, 'r') as infile, opener(output_fastq, 'w') as outfile:

        while True:
            fastq_read = list(islice(infile, 4))

            if not fastq_read:
                break

            assert len(fastq_read) == 4, 'fastq file format error'

            if not fastq_read[0].startswith('@'):
                raise ValueError('Expected @ as first character of read name')

            if not fastq_read[2].startswith('+'):
                raise ValueError('Expected @ as first character of read name')

            if tag:
                fastq_read[2] = '+' + tag + "\n"

            for fastq_line in fastq_read:
                outfile.write(fastq_line)


def merge_cell_fastqs(input_fastqs_1, output_fastq_1):
    """
    concatenates all fastq files
    if inputs are gzip: dont need to use gzip to uncompress
    and compress again since gzip files can be concatenated
    :param input_fastqs_1: input fastq file
    :type input_fastqs_1: dict or list
    :param output_fastq_1: merged fastq file
    :type output_fastq_1: str
    """
    with open(output_fastq_1, 'w') as outfile:
        for cell_id, filepath in input_fastqs_1.iteritems():
            with open(filepath) as infile:
                for line in infile:
                    outfile.write(line)


def random_subset(iterator, K):
    """ Resevoir sampling.

    Given an iterator, sample K items uniformly until
    the iterator is finished.
    """
    result = []
    N = 0

    for item in iterator:
        N += 1
        if len(result) < K:
            result.append(item)
        else:
            s = int(random.random() * N)
            if s < K:
                result[s] = item

    return result


def resample_fastqs(input_fastqs_1, input_fastqs_2, output_fastq_1, output_fastq_2, num_reads):
    """ Resample from a set of fastqs.

    Resample a fixed number of reads uniformly from a set of per cell
    paired fastq files.

    Args:
        input_fastqs_1 (dict): fastq filenames keyed by cell_id
        input_fastqs_2 (dict): fastq filenames keyed by cell_id
        output_fastq_1 (str): output fastq filename
        output_fastq_2 (str): output fastq filename
        num_reads (int): number of reads to sample
    """
    fastq_iterator = read_paired_cell_fastqs(input_fastqs_1, input_fastqs_2)
    fastq_sample = random_subset(fastq_iterator, num_reads)
    with gzip.open(output_fastq_1, 'w') as file_1, gzip.open(output_fastq_2, 'w') as file_2:
        for cell_id, fastq_lines in fastq_sample:
            for line in fastq_lines[0]:
                file_1.write(line + '\n')
            for line in fastq_lines[1]:
                file_2.write(line + '\n')


def merge_stats(input_stats, output_stats):
    """ Merge a set of destruct stats files into one file.

    Args:
        input_stats (dict): stats filenames keyed by cell_id
        output_stats (str): output stats filename
    """
    data = []
    for idx, key in enumerate(input_stats.keys()):
        data.append(pd.read_csv(input_stats[key], sep='\t'))
    data = pd.concat(data, ignore_index=True)
    data = data.groupby(['type', 'key'])['value'].sum().reset_index()
    data.to_csv(output_stats, sep='\t', index=False)


def extract_cell_counts(breakpoint_reads_filename, cell_counts_filename):
    """ Extract cell counts from a destruct breakpoint reads file

    Assumes the cell_id has been stored in the comment of the fastq and
    can be extracted from the comment column of the breakpoint reads
    table output by destruct.

    Args:
        breakpoint_reads_filename: raw breakpoint reads table from destruct
        cell_counts_filename: csv of cluster_id, cell_id, read_count
    """
    column_names = ['cluster_id', 'lib_id', 'read_id', 'read_end', 'sequence', 'quality', 'comment']
    breakpoint_reads_iter = pd.read_csv(
        breakpoint_reads_filename, sep='\t',
        names=column_names, header=None,
        iterator=True, chunksize=100000)

    cell_counts = []
    for breakpoint_reads in breakpoint_reads_iter:
        # Only count end 1 entries to remove duplicates
        breakpoint_reads = breakpoint_reads[breakpoint_reads['read_end'] == 1]
        # Comment is cell id prefixed with +
        assert breakpoint_reads['comment'].str.startswith('+').all()
        breakpoint_reads['cell_id'] = breakpoint_reads['comment'].str.lstrip('+')

        count_data = breakpoint_reads.groupby(['cluster_id', 'cell_id']).size().rename('read_count').reset_index()
        cell_counts.append(count_data)

    cell_counts = pd.concat(cell_counts, ignore_index=True)
    cell_counts.to_csv(cell_counts_filename, index=False)
