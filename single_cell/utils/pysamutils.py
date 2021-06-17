'''
Created on Jun 1, 2018

@author: dgrewal
'''
import shutil
from collections import OrderedDict

import pysam


from single_cell.utils.bamutils import bam_index


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


def _fraction_softclipped(x):
    total_softclipped = 0
    for a in x.cigar:
        if a[0] == 4:
            total_softclipped += a[1]
    return float(total_softclipped) / x.query_length


def remove_softclipped_reads(infile, outfile, softclipped_reads_threshold):
    if softclipped_reads_threshold == 1:
        shutil.copyfile(infile, outfile)
        return

    bamfile = pysam.AlignmentFile(infile, "rb")

    filteredbam = pysam.AlignmentFile(outfile, "wb", template=bamfile)
    for read in bamfile.fetch():
        if _fraction_softclipped(read) < softclipped_reads_threshold:
            filteredbam.write(read)
    filteredbam.close()

    bam_index(outfile, outfile+'.bai')
