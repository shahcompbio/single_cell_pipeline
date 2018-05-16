'''
Created on Sep 29, 2017

@author: svatrt

Updated Oct 16 2017 by dgrewal

'''
import csv
import os
import argparse
import warnings

import pandas as pd


class ConvertCSVToSEG(object):

    def __init__(
            self, segs, bin_size, metrics, output_seg, mad_threshold, multiplier):
        self.segs = segs
        self.output_seg = output_seg
        self.bin_size = bin_size
        self.metrics = metrics
        self.mad_threshold = mad_threshold
        self.multiplier = multiplier

    def get_file_format(self, infile):

        if infile.endswith('.tmp'):
            infile = infile[:-4]
        _, ext = os.path.splitext(infile)

        if ext == ".csv":
            return "csv"
        elif ext == ".gz":
            return "gzip"
        elif ext == ".h5" or ext == ".hdf5":
            return "h5"
        else:
            warnings.warn(
                "Couldn't detect output format. extension {}".format(ext))
            return "csv"

    def check_empty_file(self, path):
        """checks if file is empty
        :param path: path to the file
        :returns bool: true if file is empty, False o.w.
        """

        if not os.path.exists(path):
            raise IOError("Input file %s missing" % path)

        if os.stat(path).st_size == 0:
            return True

        with open(path) as infile:
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
        mad_cell_map = {}

        fileformat = self.get_file_format(self.metrics)

        if fileformat == "h5":
            with pd.HDFStore(self.metrics) as metrics:
                for tableid in metrics.keys():
                    multiplier = tableid.split('/')[-1]

                    if not multiplier == self.multiplier:
                        continue

                    data = metrics[tableid]
                    cellid = data["cell_id"].iloc[0]
                    mad_threshold = data["mad_neutral_state"].iloc[0]

                    assert cellid not in mad_cell_map
                    mad_cell_map[cellid] = mad_threshold

        else:
            metrics = pd.read_csv(self.metrics)
            metrics = metrics[metrics["multiplier"] == self.multiplier]

            mad_cell_map = {
                cell: mad for cell,
                mad in zip(
                    metrics["cell_id"],
                    metrics["mad_neutral_state"])}
        return mad_cell_map

    def parse_segs(self, segs, metrics):

        fileformat = self.get_file_format(segs)

        if fileformat == "h5":
            data = self.parse_segs_h5(segs, metrics)
        else:
            data = self.parse_segs_csv(segs, metrics)
        return data

    def parse_segs_h5(self, segs, metrics):

        with pd.HDFStore(segs) as segs_store:
            for tableid in segs_store.keys():

                multiplier = tableid.split('/')[-1]

                if not multiplier == self.multiplier:
                    continue

                data = segs_store[tableid]
                for _, row in data.iterrows():
                    chrom = row["chr"]
                    start = row["start"]
                    end = row["end"]
                    cell_id = row["cell_id"]
                    state = row["state"]
                    segment_length = int(end) - int(start) + 1

                    if metrics[cell_id] > self.mad_threshold:
                        continue
                    yield [cell_id, chrom, start, end, segment_length, state]

    def parse_segs_csv(self, segs, metrics):
        """parses hmmcopy segments data
        :param segs: path to hmmcopy segs file
        """

        with open(self.segs, 'r') as segfile:
            segs = csv.reader(segfile)

            lines = enumerate(segs)

            # read the header, build a dictionary with indices of headers
            _, header = next(lines)
            header = {v: i for i, v in enumerate(header)}

            # Read the segs file and write to the output file
            for _, row in lines:
                if not row["multiplier"] == self.multiplier:
                    continue

                chrom = row[header["chr"]]
                start = row[header["start"]]
                end = row[header["end"]]
                cell_id = row[header["cell_id"]]
                state = row[header["state"]]
                segment_length = int(end) - int(start) + 1

                if metrics[cell_id] > self.mad_threshold:
                    continue

                yield [cell_id, chrom, start, end, segment_length, state]

    def write_igv_segs(self, segdata, bin_width):
        """ writes IGV segs file
        :param segdata: parsed segments data,
                        format: [id, chr, start, end, width, state]
        :param bin_width: (int) extracted from reads data
        """

        with open(self.output_seg, "w") as outfile:
            # Write the first line of the output file
            outfile.write(
                '\'ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n')

            for dataval in segdata:
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

        metrics = self.read_metrics()

        segdata = self.parse_segs(self.segs, metrics)

        self.write_igv_segs(segdata, self.bin_size)


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

    parser.add_argument('--mad_threshold',
                        default=0.2,
                        type=float,
                        help='''Path to output IGV segs file''')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()

    converter = ConvertCSVToSEG(args.segs, args.bin_size,
                                args.metrics, args.output, args.mad_threshold)
    converter.main()
