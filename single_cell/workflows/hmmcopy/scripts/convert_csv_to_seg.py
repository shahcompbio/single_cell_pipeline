'''
Created on Sep 29, 2017

@author: svatrt

Updated Oct 16 2017 by dgrewal

'''
import csv
import os
import argparse


class ConvertCSVToSEG(object):

    def __init__(self, filtered_segs, filtered_reads, output_seg):
        self.filtered_segs = filtered_segs
        self.filtered_reads = filtered_reads
        self.output_seg = output_seg

    def check_empty_file(self, path):
        """checks if file is empty
        :param path: path to the file
        :returns bool: true if file is empty, False o.w.
        """

        if not os.path.exists(path):
            raise IOError("Input file %s missing" % path)

        if os.stat(self.filtered_reads).st_size != 0:
            return True

        return False

    def touch_file(self, path):
        """equivalent of posix touch to create a new empty file
        :param path: file path
        """
        open(path, "w").close()

    def get_bin_width(self, reads):
        """get the bin size from the reads data
        :param reads: path to hmmcopy reads csv file
        :returns bin_with: bin size for the reads data
        :rtype int
        """

        # read header, get index of width column, get value from next line
        with open(self.filtered_reads) as output:
            header = output.readline()
            header = header.strip().split(',')
            idx_width = header.index('width')

            line = output.readline().strip().split(',')
            bin_width = int(line[idx_width])

        return bin_width

    def parse_segs(self, segs):
        """parses hmmcopy segments data
        :param segs: path to hmmcopy segs file
        """

        with open(self.filtered_segs, 'r') as segfile:
            segs = csv.reader(segfile)

            lines = enumerate(segs)

            # read the header, build a dictionary with indices of headers
            _, header = next(lines)
            header = {v: i for i, v in enumerate(header)}

            # Read the segs file and write to the output file
            for _, row in lines:
                chrom = row[header["chr"]]
                start = row[header["start"]]
                end = row[header["end"]]
                cell_id = row[header["cell_id"]]
                state = row[header["state"]]
                segment_length = int(end) - int(start) + 1
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

                outstr = '\t'.join(dataval) + '\n'
                outfile.write(outstr)

    def main(self):

        # if the inputs are empty, then create empty output file
        if self.check_empty_file(self.filtered_reads) or\
                self.check_empty_file(self.filtered_segs):
            self.touch_file(self.output_seg)

        bin_width = self.get_bin_width(self.filtered_reads)

        segdata = self.parse_segs(self.filtered_segs)

        self.write_igv_segs(segdata, bin_width)


def parse_args():
    #=========================================================================
    # Read Command Line Input
    #=========================================================================
    parser = argparse.ArgumentParser()

    parser.add_argument('--filtered_reads',
                        required=True,
                        help='''Path to HMMcopy corrected reads output .csv file.''')

    parser.add_argument('--filtered_segs',
                        required=True,
                        help='''Path to HMMcopy segments output .csv file.''')

    parser.add_argument('--output',
                        required=True,
                        help='''Path to output IGV segs file''')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()

    converter = ConvertCSVToSEG(args.filtered_reads, args.filtered_segs,
                                args.output)
    converter.main()
