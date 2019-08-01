'''
Created on Sep 29, 2017

@author: svatrt

Updated Oct 16 2017 by dgrewal

'''
import csv
import os
import argparse
import logging
import gzip
from single_cell.utils import helpers
from single_cell.utils import csvutils

import pandas as pd
import yaml

class ConvertCSVToSEG(object):

    def __init__(
            self, segs, bin_size, metrics, output_seg, quality_threshold):
        self.segs = segs
        self.output_seg = output_seg
        self.bin_size = bin_size
        self.metrics = metrics
        self.quality_threshold = quality_threshold

    def check_empty_file(self, path):
        """checks if file is empty
        :param path: path to the file
        :returns bool: true if file is empty, False o.w.
        """

        if not os.path.exists(path):
            raise IOError("Input file %s missing" % path)

        if os.stat(path).st_size == 0:
            return True

        with helpers.getFileHandle(path) as infile:
            # header line
            _ = infile.readline()
            # data?
            data = infile.readline()
            if not data:
                return True

        return False

    def touch_file(self, path):
        """equivalent of posix touch to create a new empty file
        :param path: file path
        """
        open(path, "w").close()

    def read_metrics(self):
        """
        read metrics and get cell to mad mapping
        """
        metrics = csvutils.read_csv_and_yaml(self.metrics)

        metrics = metrics.set_index("cell_id")
        cell_order = metrics.order.sort_values().index

        # assume all cells are good, dont filter
        if 'quality' not in metrics.columns.values:
            logging.getLogger("single_cell.hmmcopy.igv_seg").warn(
                "quality column missing in data"
            )
            metrics['quality'] = 1

        qual_cell_map = {
            cell: mad for cell,
            mad in zip(
                metrics.index,
                metrics["quality"])}
        return qual_cell_map, cell_order

    def parse_segs(self, segs, metrics):
        """parses hmmcopy segments data
        :param segs: path to hmmcopy segs file
        """
        header_flag, dtypes, columns = csvutils.get_metadata(segs)

        header = {v: i for i, v in enumerate(columns)}

        segs_data = {}

        with helpers.getFileHandle(segs) as segfile:

            if header_flag:
                assert segfile.readline().strip().split(',') == columns

            for row in segfile:
                row = row.strip().split(',')

                chrom = row[header["chr"]]
                start = row[header["start"]]
                end = row[header["end"]]
                cell_id = row[header["cell_id"]]
                state = row[header["state"]]
                # float to handle scientific notation
                segment_length = int(float(end)) - int(float(start)) + 1

                if metrics[cell_id] > self.quality_threshold:
                    continue

                segs_data[cell_id] = [cell_id, chrom, start, end, segment_length, state]
        return segs_data


    def write_igv_segs(self, segdata, bin_width, cell_order):
        """ writes IGV segs file
        :param segdata: parsed segments data,
                        format: [id, chr, start, end, width, state]
        :param bin_width: (int) extracted from reads data
        """

        with open(self.output_seg, "w") as outfile:
            # Write the first line of the output file
            outfile.write(
                '\'ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n')

            for cell in cell_order:
                if cell not in segdata:
                    continue
                cell_data = segdata[cell]

                for dataval in cell_data:
                    dataval[4] = str(dataval[4] / bin_width)

                    dataval = map(str, dataval)

                    outstr = '\t'.join(dataval) + '\n'
                    outfile.write(outstr)

    def write_header(self, path):
        header = '\'ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n'

        with open(path, 'w') as outfile:
            outfile.write(header)

    def main(self):

        # if the inputs are empty, then create empty output file
        if self.check_empty_file(self.segs):
            self.write_header(self.output_seg)
            return

        qual_metrics, cell_order = self.read_metrics()

        segdata = self.parse_segs(self.segs, qual_metrics)

        self.write_igv_segs(segdata, self.bin_size, cell_order)


def parse_args():
    #=========================================================================
    # Read Command Line Input
    #=========================================================================
    parser = argparse.ArgumentParser()

    parser.add_argument('--bin_size',
                        required=True,
                        type=int,
                        help='''hmmcopy binsize''')

    parser.add_argument('--segs',
                        required=True,
                        help='''Path to HMMcopy segments output .csv file.''')

    parser.add_argument('--metrics',
                        required=True,
                        help='''Path to HMMcopy segments output .csv file.''')

    parser.add_argument('--output',
                        required=True,
                        help='''Path to output IGV segs file''')

    parser.add_argument('--quality_threshold',
                        default=0.2,
                        type=float,
                        help='''Path to output IGV segs file''')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()

    converter = ConvertCSVToSEG(args.segs, args.bin_size,
                                args.metrics, args.output, args.quality_threshold)
    converter.main()
