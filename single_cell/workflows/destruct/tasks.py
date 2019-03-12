import shutil
import random
import gzip
import pandas as pd


def merge_fastqs(input_fastqs_1, input_fastqs_2, output_fastq_1, output_fastq_2):
    """ Merge a set of pairs of fastqs into a single pair.
    """
    assert set(input_fastqs_1.keys()) == set(input_fastqs_2.keys())
    for idx, key in enumerate(input_fastqs_1.keys()):
        if idx == 0:
            shutil.copyfile(input_fastqs_1[key], output_fastq_1)
            shutil.copyfile(input_fastqs_2[key], output_fastq_2)
        else:
            with open(input_fastqs_1[key], 'rb') as infile, open(output_fastq_1, 'ab') as outfile:
                shutil.copyfileobj(infile, outfile)
            with open(input_fastqs_2[key], 'rb') as infile, open(output_fastq_2, 'ab') as outfile:
                shutil.copyfileobj(infile, outfile)


def random_subset(iterator, K):
    """ Resevoir sampling.
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


def read_paired_fastqs(input_fastqs_1, input_fastqs_2):
    """ Generator for fastq entries from paired fastq files.
    """
    assert set(input_fastqs_1.keys()) == set(input_fastqs_2.keys())
    for idx, key in enumerate(input_fastqs_1.keys()):
        opener_1 = (open, gzip.open)[input_fastqs_1[key].endswith('.gz')]
        opener_2 = (open, gzip.open)[input_fastqs_2[key].endswith('.gz')]
        with opener_1(input_fastqs_1[key], 'r') as file_1, opener_2(input_fastqs_2[key], 'r') as file_2:
            fastq_lines = [[], []]
            for fastq_1_line, fastq_2_line in zip(file_1, file_2):
                fastq_lines[0].append(fastq_1_line)
                fastq_lines[1].append(fastq_2_line)
                if len(fastq_lines[0]) == 4:
                    yield fastq_lines
                    fastq_lines = [[], []]
    assert len(fastq_lines[0]) == 0 or len(fastq_lines[0]) == 4
    if len(fastq_lines[0]) == 4:
        yield fastq_lines


def resample_fastqs(input_fastqs_1, input_fastqs_2, output_fastq_1, output_fastq_2, num_reads):
    """ Resample from a set of fastqs.
    """
    fastq_iterator = read_paired_fastqs(input_fastqs_1, input_fastqs_2)
    fastq_sample = random_subset(fastq_iterator, num_reads)
    with gzip.open(output_fastq_1, 'w') as file_1, gzip.open(output_fastq_2, 'w') as file_2:
        for fastq_lines in fastq_sample:
            for line in fastq_lines[0]:
                file_1.write(line)
            for line in fastq_lines[1]:
                file_2.write(line)


def merge_stats(input_stats, output_stats):
    """ Merge a set of destruct stats files into one file.
    """
    data = []
    for idx, key in enumerate(input_stats.keys()):
        data.append(pd.read_csv(input_stats[key], sep='\t'))
    data = pd.concat(data, ignore_index=True)
    data = data.groupby(['type', 'key'])['value'].sum().reset_index()
    data.to_csv(output_stats, sep='\t', index=False)

            

