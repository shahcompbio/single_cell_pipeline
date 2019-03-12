import shutil
import random
import gzip
import pandas as pd


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
        with opener_1(input_fastqs_1[cell_id], 'r') as file_1, opener_2(input_fastqs_2[cell_id], 'r') as file_2:
            fastq_lines = [[], []]
            for fastq_1_line, fastq_2_line in zip(file_1, file_2):
                fastq_lines[0].append(fastq_1_line.rstrip())
                fastq_lines[1].append(fastq_2_line.rstrip())
                if len(fastq_lines[0]) == 4:
                    yield fastq_lines
                    fastq_lines = [[], []]
    assert len(fastq_lines[0]) == 0 or len(fastq_lines[0]) == 4
    if len(fastq_lines[0]) == 4:
        yield cell_id, fastq_lines


def merge_tag_cell_fastqs(input_fastqs_1, input_fastqs_2, output_fastq_1, output_fastq_2):
    """ Merge a set of pairs of fastqs into a single pair.

    Merge per cell paired fastq files and modify each pair
    of reads to have a name that starts with the cell id

    Args:
        input_fastqs_1 (dict): fastq filenames keyed by cell_id
        input_fastqs_2 (dict): fastq filenames keyed by cell_id
        output_fastq_1 (str): output fastq filename
        output_fastq_2 (str): output fastq filename
    """
    fastq_iterator = read_paired_cell_fastqs(input_fastqs_1, input_fastqs_2)
    with gzip.open(output_fastq_1, 'w') as file_1, gzip.open(output_fastq_2, 'w') as file_2:
        for cell_id, fastq_lines in fastq_iterator:
            for read_end in (0, 1):
                if fastq_lines[read_end][0][0] != '@':
                    raise ValueError('Expected @ as first character of read name')
                fastq_lines[read_end][0] = '@' + cell_id + '_' + fastq_lines[read_end][0][1:]
            for line in fastq_lines[0]:
                file_1.write(line + '\n')
            for line in fastq_lines[1]:
                file_2.write(line + '\n')


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


