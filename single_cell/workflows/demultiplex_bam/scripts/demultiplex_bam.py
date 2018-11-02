import os
import pysam
import argparse
import warnings
import pandas as pd
import logging

from logging.config import fileConfig


class DemultiplexBam(object):
    """
    Demultiplexes bam files by CB (cell identifier tag) tag.

    """

    def __init__(self, ibam, obams, barcode_csv):
        """

        :param ibam: bam to demultiplex
        :param obams: dictionary of output bam file names
        :param barcode_csv: csv file containing barcodes of aligned reads to demultiplex.
        barcodes must be in the fist column of the csv file.
        """

        self.bam = self.__get_bam_reader(ibam)

        self.obams = obams

        self.barcodes = self.__get_barcodes(barcode_csv)

    def __get_barcodes(self, barcode_file):

        df = pd.read_csv(barcode_file)
        return set(df.barcode.tolist())

    def __enter__(self):

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):

        self.bam.close()

    def __get_bam_reader(self, bam_file):
        """returns pysam bam object
        :returns pysam bam object
        """

        return pysam.AlignmentFile(bam_file, 'rb')

    def __fetch(self, chrom, start, end):
        """returns iterator over reads in the specified region
        :param start: bin starting pos (int)
        :param end: bin end pos (int)
        :returns iterator over reads
        """

        return self.bam.fetch(chrom, start, end)

    def __get_file_for_tag(self, tag_to_file_map, tag_value):
        """
        Get a file handle for the specified tag_value.

        :param tag_to_file_map: map from tag value to file handler that is used for caching file handles
        :param tag_value: tag value
        :return: File handle for the specified tag value
        """

        if tag_value not in tag_to_file_map.keys():
            file_handle = pysam.AlignmentFile(tag_value, "wb", template=self.bam)
            tag_to_file_map[tag_value] = file_handle

        return tag_to_file_map[tag_value]

    def demultiplex_bam(self):
        """
        Demultiplexes the input bam.

        """

        read_count = 0
        tag_to_file_map = {}

        for read in self.bam.fetch():
            read_count += 1
            tags = dict(read.tags)
            if 'CB' in tags.keys():
                if tags['CB'] in self.barcodes:  # filter by barcodes
                    outfile = self.__get_file_for_tag(tag_to_file_map, self.obams[tags['CB']])
                    outfile.write(read)
                else:
                    pass  # discard read
            else:
                outfile = self.__get_file_for_tag(tag_to_file_map, self.obams["undetermined"])
                outfile.write(read)

        logging.info("closing demultiplexed bam files")
        for file in tag_to_file_map.values():
            try:
                file.close()
            except:
                logging.error('error closing file ' + file)

        logging.info("num reads demultiplexed = " + str(read_count))



