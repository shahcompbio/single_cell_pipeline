'''
Created on Nov 20, 2015

@author: Andrew Roth
'''
from collections import namedtuple

import csv
import pysam

from components_utils import flatten_input

chrom_map = {'X': 23, 'Y': 24, 'M': 25, 'MT': 25}


def merge_vcfs(in_files, out_file):
    in_files = flatten_input(in_files)

    with open(out_file, 'w') as out_fh:
        write_header(out_fh)

        writer = csv.DictWriter(out_fh, ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'], delimiter='\t')

        reader = MultiVcfReader(in_files)

        for row in reader:
            writer.writerow({'CHROM': row.chrom,
                             'POS': row.coord,
                             'ID': '.',
                             'REF': row.ref,
                             'ALT': row.alt,
                             'QUAL': '.',
                             'FILTER': '.',
                             'INFO': '.'})

        reader.close()


def get_chrom_order(chrom):
    '''
    Convert chromosome names so they will sort 1, 2, 3, ..., X, Y, MT, etc..
    '''
    chrom = chrom.replace('chr', '')

    chrom = chrom.replace('Chr', '')

    try:
        chrom = int(chrom)

    except ValueError:
        if chrom in chrom_map:
            chrom = chrom_map[chrom]

        else:
            chrom = chrom

    return chrom


def write_header(fh):
    fh.write('##fileformat=VCFv4.1\n')

    header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    header = '\t'.join(header)

    fh.write('#{0}\n'.format(header))

LightVCFRecord = namedtuple('LightVCFRecord', ['chrom', 'coord', 'ref', 'alt'])


class MultiVcfReader(object):

    def __init__(self, vcf_files):
        self._readers = []

        for file_name in vcf_files:
            self._readers.append(pysam.Tabixfile(file_name, parser=pysam.asVCF()))

    def __iter__(self):
        for chrom in self.chroms:
            active_iters = self._load_iters(chrom)

            while len(active_iters) > 0:
                pos_buffer = set()

                min_coord = min([x.record.coord for x in active_iters])

                new_active_iters = []

                for chrom_iter in active_iters:
                    while chrom_iter.coord == min_coord:
                        record = chrom_iter.record

                        # Handles multiple alt alleles.
                        for alt in record.alt.split(','):
                            pos_buffer.add(LightVCFRecord(record.chrom, record.coord, record.ref, alt))

                        chrom_iter.buffer_next_record()

                    if chrom_iter.active:
                        new_active_iters.append(chrom_iter)

                for record in sorted(pos_buffer, key=lambda x: (get_chrom_order(x.chrom), x.coord, x.ref, x.alt)):
                    yield record

                active_iters = new_active_iters

    def close(self):
        for reader in self._readers:
            reader.close()

    @property
    def chroms(self):
        '''
        Get a union set of chromosomes present in VCF readers.
        '''
        chroms = set()

        for reader in self._readers:
            chroms.update(set(reader.contigs))

        return sorted(chroms, key=lambda x: get_chrom_order(x))

    def _load_iters(self, chrom):
        iters = []

        for reader in self._readers:
            try:
                chrom_iter = reader.fetch(chrom)

            except (KeyError, ValueError):
                continue

            iters.append(BufferedChromIter(chrom_iter))

        return iters


class BufferedChromIter(object):

    def __init__(self, chrom_iter):
        self._iter = chrom_iter

        self.active = True

        self.buffer_next_record()

    def buffer_next_record(self):
        try:
            record = self._iter.next()

            self.record = LightVCFRecord(record.contig, record.pos + 1, record.ref, record.alt)

            self.coord = self.record.coord

        except StopIteration:
            self.active = False

            self.record = None

            self.coord = None
