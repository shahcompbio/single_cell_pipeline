'''
Created on Jul 14, 2017

@author: dgrewal
'''
from collections import Counter

import csv
import pysam
import gzip

def main(args):
    sampid = args.sample_id
    
    bam_file = pysam.Samfile(args.bam_file, 'rb')

    bases, fields = get_bases_and_fields(args)

    writer = csv.DictWriter(open(args.out_file, 'w'), fields, delimiter=',')

    writer.writeheader()

    positions, target_bases = load_positions(args.positions_file)

    pileup_iterator = positions_iterator(
        bam_file, positions, args.max_pileup_depth)

    for pileup_column in pileup_iterator:
        counts = get_counts(pileup_column,
                            args.min_bqual,
                            args.min_mqual,
                            args.count_duplicates,
                            args.count_qc_failures,
                            args.strand_counts)

        total_counts = sum([counts[b] for b in bases])

        if total_counts < args.min_counts:
            continue

        out_row = {}

        out_row['chrom'] = bam_file.getrname(pileup_column.tid)

        # One based coordinate
        out_row['coord'] = pileup_column.pos + 1

        pos = (out_row['chrom'], out_row['coord'])

        ref_base = target_bases[pos]['ref_base']

        var_base = target_bases[pos]['var_base']

        out_row['ref_base'] = ref_base

        out_row['var_base'] = var_base

        if args.strand_counts:
            out_row[sampid+'_ref_counts_forward'] = counts[ref_base.upper()]

            out_row[sampid+'_ref_counts_reverse'] = counts[ref_base.lower()]

            out_row[sampid+'_var_counts_forward'] = counts[var_base.upper()]

            out_row[sampid+'_var_counts_reverse'] = counts[var_base.lower()]

        else:
            out_row[sampid+'_ref_counts'] = counts[ref_base]

            out_row[sampid+'_var_counts'] = counts[var_base]

        writer.writerow(out_row)


def get_bases_and_fields(args):
    sampid = args.sample_id

    if args.strand_counts:
        bases = ['A', 'a', 'C', 'c', 'G', 'g', 'T', 't']
    else:
        bases = ['A', 'C', 'G', 'T']


    fields = ['chrom', 'coord', 'ref_base', 'var_base']

    if args.strand_counts:
        fields += [sampid+'_ref_counts_forward', sampid+'_ref_counts_reverse',
                   sampid+'_var_counts_forward', sampid+'_var_counts_reverse']

    else:
        fields += [sampid+'_ref_counts', sampid+'_var_counts']

    return bases, fields


def get_counts(pileup_column, min_bqual, min_mqual, count_duplicates, count_qc_failures, strand_counts):
    bases = []

    for pileup_read in pileup_column.pileups:
        # Skip flagged duplicates if we don't want them
        if pileup_read.alignment.is_duplicate and not count_duplicates:
            continue

        # Skip QC failures if we don't want them
        if pileup_read.alignment.is_qcfail and not count_qc_failures:
            continue

        if pileup_read.is_del:
            continue

        mqual = pileup_read.alignment.mapq

        # Skip positions with low mapping quality
        if mqual < min_mqual:
            continue

        # Nucleotide handling
        else:
            bqual = ord(pileup_read.alignment.qual[pileup_read.query_position]) - 33

            if bqual < min_bqual:
                continue

            base = pileup_read.alignment.seq[pileup_read.query_position]

            if strand_counts:
                if pileup_read.alignment.is_reverse:
                    base = base.lower()

                else:
                    base = base.upper()

            else:
                base = base.upper()

            bases.append(base)

    return Counter(bases)


def is_gzip(filename):
    """
    Uses the file contents to check if the file is gzip or not.
    The magic number for gzip is 1f 8b
    See KRONOS-8 for details
    """
    with open(filename) as f:
        file_start = f.read(4)
    
        if file_start.startswith("\x1f\x8b\x08"):
            return True
        return False

def get_var_base(counts, ref_base):
    # Get rid of strand information
    counts = Counter([x.upper() for x in counts.elements()])

    # Remove reference base
    del counts[ref_base]

    if len(counts) == 0:
        var_base = 'N'

    else:
        # Get most common non-reference base
        var_base, var_counts = counts.most_common()[0]

    return var_base


def load_positions(file_name):

    positions = []
    target_bases = {}


    if is_gzip(file_name):
        merge_file = gzip.open(file_name, 'rb')
    else:
        merge_file = open(file_name)

    header = merge_file.readline().strip().split()
    chr_idx = header.index('chromosome')
    pos_idx = header.index('start')
    ref_idx = header.index('ref')
    alt_idx = header.index('alt')



    for line in merge_file:
        columns = line.split("\t")
        chrom = columns[chr_idx]
        pos = int(columns[pos_idx])
        ref = columns[ref_idx]
        alt  = columns[alt_idx]

        pos = (chrom, pos)
        target_bases[pos] = {'ref_base': ref, 'var_base': alt}
        
        positions.append(pos)
    
    merge_file.close()

    return positions, target_bases



def positions_iterator(bam_file, positions, max_pileup_depth):
    for chrom, coord in positions:
        if bam_file.count(chrom, coord - 1, coord) == 0:
            pos = coord - 1

            tid = bam_file.gettid(chrom)

            pileup_column = DummyPileupProxy(pos, tid)

        else:
            pileup_iterator = bam_file.pileup(chrom,
                                              coord - 1,
                                              coord,
                                              truncate=True,
                                              max_depth=max_pileup_depth,
                                              mask=False,
                                              stepper='all')

            for pileup_column in pileup_iterator:
                break

        yield pileup_column


class DummyPileupProxy(object):

    def __init__(self, pos, tid):
        self.pileups = []

        self.pos = pos

        self.tid = tid

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('bam_file',
                        help='''BAM file to generate counts from.''')

    parser.add_argument('positions_file',
                        help='''vcf file which lists the positions to get counts from. File must
                        contain `chrom` and `coord` fields which specify the chromsome and chromosome coordinates
                        respectively for the positions''')

    parser.add_argument('out_file',
                        help='''Path where counts file will be written.''')

    parser.add_argument('sample_id',
                        help='''identifier for the sample''')

    parser.add_argument('--count_duplicates', default=False, action='store_true',
                        help='''If set flagged duplicates will be counted. Default is to ignore duplicates.''')

    parser.add_argument('--count_qc_failures', default=False, action='store_true',
                        help='''If set flagged QC failures will be counted. Default is to ignore QC failures.''')

    parser.add_argument('--max_pileup_depth', type=int, default=int(1e7),
                        help='''Max coverage depth of a position pileup engine will consider. Default is 1e7.''')

    parser.add_argument('--min_bqual', type=int, default=0,
                        help='''Minimum base quality required to count a read. Default 0.''')

    parser.add_argument('--min_counts', type=int, default=0,
                        help='''Minimum number of counts required to report a site. Default is 0.''')

    parser.add_argument('--min_mqual', type=int, default=0,
                        help='''Minimum mapping quality required to count a read. Default 0.''')

    parser.add_argument('--strand_counts', default=False, action='store_true',
                        help='''If set the strand of each base will be counted. Instead of four fields (A, C, G, T),
                        eight fields (A, a, C, c, G, g, T, t) will be reported where upper case values are the forward
                        strand and lower case values are the reverse strand. Default is to ignore strand.''')

    args = parser.parse_args()

    main(args)