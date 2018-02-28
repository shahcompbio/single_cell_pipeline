'''
Created on Oct 31, 2015

@author: Andrew Roth
'''
from collections import OrderedDict

import pysam
import vcf
import pandas as pd


default_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']


def get_regions(chromosome_lengths, split_size):
    if split_size is None:
        return dict(enumerate(chromosome_lengths.keys()))

    regions = {}
    region_index = 0

    for chrom, length in chromosome_lengths.iteritems():
        lside_interval = range(1, length + 1, split_size)
        rside_interval = range(split_size, length + split_size, split_size)

        for beg, end in zip(lside_interval, rside_interval):
            end = min(end, length)

            regions[region_index] = '{}:{}-{}'.format(chrom, beg, end)
            region_index += 1

    return regions


def get_vcf_regions(vcf_file, split_size, chromosomes=None):
    if split_size is None:
        return dict(enumerate(chromosomes))
    chromosome_lengths = load_vcf_chromosome_lengths(vcf_file, chromosomes=chromosomes)
    return get_regions(chromosome_lengths, split_size)


def get_bam_regions(bam_file, split_size, chromosomes=None):
    chromosome_lengths = load_bam_chromosome_lengths(bam_file, chromosomes=chromosomes)
    return get_regions(chromosome_lengths, split_size)


def calculate_vcf_chromosome_lengths(file_name, chromosomes=None):
    if file_name.endswith('gz'):
        compression = 'gzip'
    else:
        compression = None

    tsv_reader = pd.read_csv(
        file_name, sep='\t', comment='#', chunksize=int(1e6),
        names=['chrom', 'coord'], usecols=[0, 1], index_col=0,
        converters={'chrom': str}, compression=compression)

    max_coord = [chunk.groupby(level=0).max() for chunk in tsv_reader]
    max_coord = pd.concat(max_coord).groupby(level=0).max()

    chromosome_lengths = max_coord + 1000

    return chromosome_lengths['coord'].to_dict()


def load_vcf_chromosome_lengths(file_name, chromosomes=None):
    chromosome_lengths = OrderedDict()

    vcf_reader = vcf.Reader(filename=file_name)

    if len(vcf_reader.contigs) == 0:
        return calculate_vcf_chromosome_lengths(file_name, chromosomes=chromosomes)

    if chromosomes is None:
        chromosomes = vcf_reader.contigs.keys()

    else:
        chromosomes = chromosomes

    calc_lens = calculate_vcf_chromosome_lengths(file_name, chromosomes=chromosomes)

    for chrom, contig in vcf_reader.contigs.iteritems():
        assert chrom == contig.id

        if chrom not in chromosomes:
            continue

        if contig.length is None:
            chromosome_lengths[str(chrom)] = calc_lens[str(chrom)]

        else:
            chromosome_lengths[str(chrom)] = int(contig.length)

    if len(chromosome_lengths) == 0:
        raise Exception('no chromosomes found in vcf header')

    return chromosome_lengths


def load_bam_chromosome_lengths(file_name, chromosomes=None):
    chromosome_lengths = OrderedDict()

    bam = pysam.Samfile(file_name, 'rb')

    if chromosomes is None:
        chromosomes = bam.references

    else:
        chromosomes = chromosomes

    for chrom, length in zip(bam.references, bam.lengths):
        if chrom not in chromosomes:
            continue

        chromosome_lengths[str(chrom)] = int(length)

    return chromosome_lengths


def parse_region_for_vcf(region):
    if ':' not in region:
        return region, None, None

    chrom, coords = region.split(':')

    if '-' not in coords:
        return chrom, int(coords) - 1, None

    beg, end = coords.split('-')

    beg = int(beg) - 1
    end = int(end)

    return chrom, beg, end
