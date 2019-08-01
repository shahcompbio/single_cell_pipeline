'''
Created on Jun 1, 2018

@author: dgrewal
'''
from collections import OrderedDict
import pysam

def load_chromosome_lengths(file_name, chromosomes=None):

    chromosome_lengths = OrderedDict()

    ref = pysam.Fastafile(file_name)

    for chrom, length in zip(ref.references, ref.lengths):
        if chromosomes and chrom not in chromosomes:
            continue

        chromosome_lengths[str(chrom)] = int(length)

    return chromosome_lengths

def get_regions_from_reference(reference_fastq, split_size, chromosomes):

    chromosome_lengths = load_chromosome_lengths(
        reference_fastq,
        chromosomes=chromosomes
    )
    return get_regions(chromosome_lengths, split_size)


def get_regions(chromosome_lengths, split_size):
    if split_size is None:
        return dict(enumerate(chromosome_lengths.keys()))

    regions = []

    for chrom, length in chromosome_lengths.items():
        lside_interval = range(1, length + 1, split_size)
        rside_interval = range(split_size, length + split_size, split_size)

        for beg, end in zip(lside_interval, rside_interval):
            end = min(end, length)

            regions.append('{}-{}-{}'.format(chrom, beg, end))

    return regions
