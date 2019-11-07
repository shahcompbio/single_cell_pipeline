import gzip
import os
import random
import shutil
from itertools import islice

import pandas as pd
from single_cell.utils import helpers

import pypeliner

from single_cell.utils import fastqutils



def destruct_bamdisc_and_numreads(
        destruct_config, normal_bam_file, stats,
        reads_1, reads_2, sample_1, sample_2, tempdir
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

    cmd = ['destruct_bamdiscordantfastq',
           '-r',
           '-c', destruct_config['bam_max_soft_clipped'],
           '-f', destruct_config['bam_max_fragment_length'],
           '-b', normal_bam_file,
           '-s', stats,
           '--fastq1', reads_1,
           '--fastq2', reads_2,
           '-t', tempdir,
           '-n', destruct_config['num_read_samples'],
           '--sample1', sample_1,
           '--sample2', sample_2,
           ]
    pypeliner.commandline.execute(*cmd)

    numreads_r1 = get_read_count(reads_1)
    numreads_r2 = get_read_count(reads_2)

    return max(numreads_r1, numreads_r2)


def merge_fastqs(inputs, output):
    read_counter = 0
    with helpers.getFileHandle(output, 'wt') as merged:
        for cellid in inputs:
            infile = inputs[cellid]
            reader = fastqutils.FastqReader(infile)
            for read in reader.get_read_iterator():
                read[0] = '@' + str(int(read_counter)) + '/' + read[0].split('/')[1]
                for line in read:
                    merged.write(line)
                read_counter += 1




def merge_read_counts(readcounts):
    return readcounts


def re_index_reads_both(
        input_r1, reindex_r1, input_r2, reindex_r2,
        cell_id, cells, offset, tag=False
):
    re_index_reads(input_r1, reindex_r1, cell_id, cells, offset, tag=tag)
    re_index_reads(input_r2, reindex_r2, cell_id, cells, offset, tag=tag)


def get_start_count(cells, readcounts, cell_id):
    index = cells.index(cell_id)

    reads_before = [readcounts[cell] for cell in cells[:index]]

    return sum(reads_before)


def re_index_reads(input_fastq, output_fastq, cell_id, cells, cell_read_counts, tag=False):
    start_count = get_start_count(cells, cell_read_counts, cell_id)

    with helpers.getFileHandle(input_fastq) as infile:
        with helpers.getFileHandle(output_fastq, 'wt') as outfile:

            while True:
                fastq_read = list(islice(infile, 4))

                if not fastq_read:
                    break

                assert len(fastq_read) == 4, 'fastq file format error'

                if not fastq_read[0].startswith('@'):
                    raise ValueError('Expected @ as first character of read name')

                if not fastq_read[2].startswith('+'):
                    raise ValueError('Expected = as first character of read comment')

                fastq_read[0] = '@' + str(int(start_count)) + '/' + fastq_read[0].split('/')[1]

                start_count += 1

                if tag:
                    fastq_read[2] = '+' + cell_id + "\n"

                for fastq_line in fastq_read:
                    outfile.write(fastq_line)


def get_read_count(input_fastq):
    def blocks(files, size=65536):
        while True:
            b = files.read(size)
            if not b:
                break
            yield b

    with helpers.getFileHandle(input_fastq) as indata:
        linecount = sum(bl.count("\n") for bl in blocks(indata))

    assert linecount % 4 == 0

    readcount = linecount / 4

    return int(readcount)


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
    with open(output_fastq_1, 'wb') as outfile:
        for cell_id, filepath in input_fastqs_1.items():
            with open(filepath, 'rb') as infile:
                shutil.copyfileobj(infile, outfile, length=16 * 1024 * 1024)


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


def read_paired_cell_fastqs(input_fastqs_1, input_fastqs_2):
    """ Generator for fastq entries from paired fastq files.
    Read per cell paired fastq files and for each entry in the
    matching fastq pair, return the cell id and a list of
    two lists of fastq lines.
    Args:
        input_fastqs_1 (dict): fastq filenames keyed by cell_id
        input_fastqs_2 (dict): fastq filenames keyed by cell_id
    Yields:
        cell_id (str), fastq_lines (list)
    """
    assert set(input_fastqs_1.keys()) == set(input_fastqs_2.keys())
    for idx, cell_id in enumerate(input_fastqs_1.keys()):
        opener_1 = (open, gzip.open)[input_fastqs_1[cell_id].endswith('.gz')]
        opener_2 = (open, gzip.open)[input_fastqs_2[cell_id].endswith('.gz')]
        with opener_1(input_fastqs_1[cell_id], 'rt') as file_1, opener_2(input_fastqs_2[cell_id], 'rt') as file_2:
            fastq_lines = [[], []]
            for fastq_1_line, fastq_2_line in zip(file_1, file_2):
                fastq_lines[0].append(fastq_1_line.rstrip())
                fastq_lines[1].append(fastq_2_line.rstrip())
                if len(fastq_lines[0]) == 4:
                    yield cell_id, fastq_lines
                    fastq_lines = [[], []]
    assert len(fastq_lines[0]) == 0 or len(fastq_lines[0]) == 4, fastq_lines
    if len(fastq_lines[0]) == 4:
        yield cell_id, fastq_lines


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
    with gzip.open(output_fastq_1, 'wt') as file_1, gzip.open(output_fastq_2, 'wt') as file_2:
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
    column_names = ['cluster_id', 'lib_id', 'read_id', 'read_end', 'sequence', 'quality', 'comment', 'filtered']
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


def filter_reads_file(
        breakpoint_reads_filename, destruct_table,
        filtered_breakpoint_reads_filename
):
    # writing to file in append mode, so make sure that it doesnt
    # exist before we start appending
    if os.path.exists(filtered_breakpoint_reads_filename):
        os.remove(filtered_breakpoint_reads_filename)

    destruct = pd.read_csv(destruct_table, sep='\t')
    passed_brks = set(destruct['prediction_id'])

    column_names = ['cluster_id', 'lib_id', 'read_id', 'read_end', 'sequence', 'quality', 'comment', 'filtered']
    breakpoint_reads_iter = pd.read_csv(
        breakpoint_reads_filename, sep='\t',
        names=column_names, header=None,
        iterator=True, chunksize=100000)

    for i, breakpoint_reads in enumerate(breakpoint_reads_iter):
        passed_reads = breakpoint_reads[breakpoint_reads['cluster_id'].isin(passed_brks)]
        passed_reads = passed_reads[passed_reads['filtered'] == False]
        passed_reads.to_csv(filtered_breakpoint_reads_filename, sep='\t', header=False, index=False, mode='a')
