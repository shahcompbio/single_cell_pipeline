import pandas as pd

default_chromosomes = [str(a) for a in range(1, 23)] + ['X', 'Y']


def read_chromosome_lengths(genome_fasta_index, chromosomes=default_chromosomes):
    fai = pd.read_csv(genome_fasta_index, sep='\t', header=None, names=['chrom', 'length', 'V3', 'V4', 'V5'])
    fai = fai.set_index('chrom')['length']
    fai = fai.reindex(chromosomes).astype(int)
    return fai.to_dict()


def get_split_regions(split_size, refgenome, chromosomes=default_chromosomes):
    genome_fasta_index = refgenome + '.fai'

    chromosome_lengths = read_chromosome_lengths(genome_fasta_index, chromosomes=chromosomes)

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
