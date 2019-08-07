'''
Created on Oct 10, 2017

@author: dgrewal
'''
import argparse
import os

import numpy as np
import pandas as pd
import pysam
import logging

class ReadCounter(object):
    """
    calculate reads per bin from the input bam file
    """

    def __init__(
            self, bam, output, window_size, chromosomes, mapq,
            seg=None, excluded=None, reference=None
    ):
        self.bam = bam

        self.output = output

        self.window_size = window_size

        if chromosomes:
            self.chromosomes = chromosomes
        else:
            self.chromosomes = self.__get_chr_names()

        self.bam = self.__get_bam_reader()
        self.chr_lengths = self.__get_chr_lengths()

        self.mapq_threshold = mapq

        self.seg = seg

        if excluded is not None:
            self.excluded = pd.read_csv(excluded, sep="\t", )
            self.excluded.columns = ["chrom", "start", "end"]
        else:
            self.excluded = None

        self.reference = reference

    def __get_bam_header(self):
        return self.bam.header

    def __get_chrom_excluded(self, chrom, chrom_length):
        # chrom_excluded = np.zeros(chrom_length, dtype=np.uint8)
        # add some padding in case the list is 1 based and chr length is 0 based
        chrom_excluded = np.zeros(chrom_length + 1, dtype=np.uint8)

        for start, end in self.excluded.loc[self.excluded['chrom'] == chrom, ['start', 'end']].values:
            start = min(start, chrom_length)
            end = min(end, chrom_length)
            chrom_excluded[start:end] = 1

        return chrom_excluded

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # clean up output if there are any exceptions
        if exc_type and os.path.exists(self.output):
            os.remove(self.output)

    def __get_chr_lengths(self):
        """ returns dict with chromosome names and lengths
        :returns dictionary with chromosome name (str) and lengths(int)
        :rtype dictionary
        """

        names = self.bam.references
        lengths = self.bam.lengths
        return {name: length for name, length in zip(names, lengths)}

    def __get_bam_reader(self):
        """returns pysam bam object
        :returns pysam bam object
        """
        return pysam.AlignmentFile(self.bam, 'rb')

    def __get_chr_names(self):
        """extracts chromosome names from the bam file
        :returns list of chromosome names
        :rtype list
        """
        return self.bam.references

    def __fetch(self, chrom, start, end):
        """returns iterator over reads in the specified region
        :param chrom: chromosome name (str)
        :param start: bin starting pos (int)
        :param end: bin end pos (int)
        :returns iterator over reads
        """
        return self.bam.fetch(chrom, start, end)

    def filter(self, pileupobj, chrom_excluded=None):
        """remove low mapping quality reads and duplicates
        :param pileupobj: pysam read object
        :returns boolean: false if the read passes filters.
        :rtype boolean
        """

        pos = pileupobj.reference_start
        if chrom_excluded is not None and chrom_excluded[pos]:
            return True

        if pileupobj.is_duplicate:
            return True

        if pileupobj.mapping_quality < self.mapq_threshold:
            return True

        try:
            fastqscreen_tags = pileupobj.get_tag('FS')
        except KeyError:
            logging.getLogger("read_counter").warn(
                "couldn't get FS tag from bam"
            )
            fastqscreen_tags = None
        if fastqscreen_tags and self.reference:
            fastqscreen_tags = fastqscreen_tags.split(',')
            fastqscreen_tags = [val.split('_') for val in fastqscreen_tags]
            fastqscreen_tags = {val[0]: val[1] for val in fastqscreen_tags}

            if int(fastqscreen_tags[self.reference]) == 0:
                return True

        return False

    def write_header(self, chrom, outfile):
        """writes headers, single header if seg format,
        one header per chromosome otherwise.
        :param chrom: chromosome name
        :param outfile: output file object.
        """
        if self.seg:
            outfile.write("data\tchr\tstart\tend\tcount\n")
        else:
            outstr = "fixedStep chrom=%s start=1 step=%s span=%s\n" \
                     % (chrom, self.window_size, self.window_size)
            outfile.write(outstr)

    def write(self, chrom, start, stop, count, outfile):
        """writes bin and counts to the output file.
        supports seg and wig formats
        :param chrom: chromosome name
        :param start: bin start
        :param stop: bin stop
        :param count: no of reads in the bin
        :param outfile: output file object.
        """
        if self.seg:
            outstr = '\t'.join(
                map(str, ['reads', chrom, start, stop, count])) + '\n'
            outfile.write(outstr)
        else:
            outfile.write(str(count) + '\n')

    def get_data(self, data, chrom, outfile):
        """iterates over reads, calculates counts and writes to output
        :param data: pysam iterator over reads
        :param chrom: str: chromosome name
        :param outfile: output file object
        """
        reflen = self.chr_lengths[chrom]

        chrom_excluded = None
        if self.excluded is not None:
            chrom_excluded = self.__get_chrom_excluded(chrom, reflen)

        count = 0
        start = 0
        end = start + self.window_size
        for pileupobj in data:
            while pileupobj.pos > end:
                self.write(chrom, start, end, count, outfile)
                count = 0
                start += self.window_size
                end = min(start + self.window_size, reflen)

            if not self.filter(pileupobj, chrom_excluded):
                count += 1

        # catch up with end of chromosome
        while True:
            self.write(chrom, start, end, count, outfile)
            count = 0
            start += self.window_size
            end = start + self.window_size
            if start > reflen:
                break
            if end > reflen:
                self.write(chrom, start, reflen, count, outfile)
                break

    def main(self):
        """for each chromosome, iterate over all reads. use starting position
        of the read to calculate read counts per bin (no double counting).
        """
        with open(self.output, 'w') as outfile:
            if self.seg:
                self.write_header(None, outfile)

            for chrom in self.chromosomes:
                if not self.seg:
                    self.write_header(chrom, outfile)

                reflen = self.chr_lengths[chrom]

                # get read iterator for the full chromosome
                # code assumes the iterator is sorted.
                data = self.__fetch(chrom, 0, reflen)

                self.get_data(data, chrom, outfile)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('bam',
                        help='specify the path to the input bam file')

    parser.add_argument('output',
                        help='specify path to the output file')
    parser.add_argument('--chromosomes',
                        nargs='*',
                        default=map(str, range(1, 23)) + ['X', 'Y'],
                        help='specify target chromosomes'
                        )
    parser.add_argument('-w', '--window_size',
                        type=int,
                        default=1000,
                        help='specify bin size')
    parser.add_argument('-m', '--mapping_quality_threshold',
                        type=int,
                        default=0,
                        help='threshold for the mapping quality, reads ' \
                             'with quality lower than threshold will be ignored')

    parser.add_argument('--seg',
                        default=False,
                        action='store_true',
                        help='write the output in seg format')

    parser.add_argument('--exclude_list',
                        default=None,
                        help='regions to skip')

    parser.add_argument('--reference', default=None)

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = parse_args()
    with ReadCounter(args.bam, args.output, args.window_size,
                     args.chromosomes, args.mapping_quality_threshold,
                     args.seg, excluded=args.exclude_list,
                     reference=args.reference) as rcount:
        rcount.main()
