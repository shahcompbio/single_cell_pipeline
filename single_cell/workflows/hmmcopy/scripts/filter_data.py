'''
Created on Nov 16, 2016

@author: dgrewal
'''
from __future__ import division
import logging
import argparse
import pandas as pd
import math
import gzip


def parse_args():
    #=========================================================================
    # Read Command Line Input
    #=========================================================================
    parser = argparse.ArgumentParser()

    parser.add_argument('--corrected_reads',
                        required=True,
                        help='''Path to HMMcopy corrected reads output .csv file.''')

    parser.add_argument('--segments',
                        required=True,
                        help='''Path to HMMcopy segments output .csv file.''')

    parser.add_argument('--quality_metrics',
                        required=True,
                        help='''Optional quality metrics file for the run, with 'mad_neutral_state' column.''')

    parser.add_argument('--mad_threshold', type=float, default=0,
                        help='''all cells that have low MAD won't be plotted''')

    parser.add_argument('--reads_output',
                        required=True,
                        help='''Path to HMMcopy corrected reads output .pdf file.''')

    parser.add_argument('--segs_output',
                        required=True,
                        help='''Path to HMMcopy segs reads output .pdf file.''')

    args = parser.parse_args()
    return args


class FilterHmmData(object):
    """
    generate the reads, bias and segment plots
    """

    def __init__(self, quality_metrics, segments, reads,
                 mad_threshold, reads_out, segments_out, compression=None):

        self.quality_metrics = quality_metrics
        self.segments = segments
        self.reads = reads
        self.mad_threshold = mad_threshold
        self.reads_out = reads_out
        self.segments_out = segments_out

        if compression and not compression=="gzip":
            raise Exception("Cannot compress as {}, only supports gzip".format(compression))

        self.compression = compression

        if not self.mad_threshold:
            self.mad_threshold = 0

    def load_data_pandas(self, infile):
        """

        """
        data = pd.read_csv(infile,
                           sep=',')

        data = data.groupby('cell_id')

        return data

    def read_quality_metrics(self):
        """

        """

        df = self.load_data_pandas(self.quality_metrics)

        return df

    def get_mad_score(self, sample_id, metrics):
        """
        """
        mad = metrics.get_group(sample_id)['mad_neutral_state'].iloc[0]
        return mad

    def check_mad_score(self, sample, metrics):
        """

        """
        mad = self.get_mad_score(sample, metrics)

        # if mad_threshold is set to nonzero.
        # zero is defaults and means mad_threshold is not set. so no filtering
        if self.mad_threshold:
            if math.isnan(mad):
                return False

            if mad > self.mad_threshold:
                return False
        return True

    def filter_write(self, infile, metrics, outfile):

        with open(infile) as inp:

            if self.compression == "gzip":
                output = gzip.open(outfile, 'w')
            else:
                output = open(outfile, 'w')

            header = inp.readline()
            #if input file is empty
            if not header:
                logging.getLogger("single_cell.hmmcopy.filter_data").warn("no data to filter")
                return

            output.write(header)
            header = header.strip().split(',')
            samp_idx = header.index('cell_id')

            for line in inp:
                samp = line.strip().split(',')[samp_idx]
                if self.check_mad_score(samp, metrics):
                    output.write(line)

            output.close()

    def main(self):
        """
        main
        """
        metrics = self.read_quality_metrics()

        self.filter_write(self.reads, metrics, self.reads_out)
        self.filter_write(self.segments, metrics, self.segments_out)

if __name__ == '__main__':
    args = parse_args()

    genhmm = FilterHmmData(args.quality_metrics, args.segments,
                           args.corrected_reads, args.mad_threshold,
                           args.reads_output, args.segs_output)

    genhmm.main()
