'''
Created on Apr 13, 2018

@author: dgrewal
'''
import numpy as np
import argparse
import warnings


def read_wig(infile, counts=False):
    """read wiggle files

    :param infile: input wiggle file
    :param counts: set to true if infile wiggle has integer values
    """

    data = {}

    with open(infile) as wig:
        for line in wig:
            line = line.strip()

            if line.startswith('fixedStep'):
                line = line.strip().split()

                chrom = line[1].split('=')[1]
                start = int(line[2].split('=')[1])
                step = int(line[3].split('=')[1])
                span = int(line[4].split('=')[1])
                chrom_key = (chrom, start, step, span)
            else:
                value = int(line) if counts else float(line)

                if chrom_key not in data:
                    data[chrom_key] = [value]
                else:
                    data[chrom_key].append(value)

    return data


def merge_two_wigs(wig1, wig2):
    merged_wig = {}

    chromosomes = sorted(set(wig1.keys()).union( set(wig2.keys())))

    for chrom in chromosomes:
        wig1data = wig1.get(chrom, None)
        wig2data = wig2.get(chrom, None)

        if not isinstance(wig2data, np.ndarray) and not wig2data:
            #just add zeros if no data
            wig2data = [0] * len(wig1data)

        if not isinstance(wig1data, np.ndarray) and not wig1data:
            #just add zeros if no data
            wig1data = [0] * len(wig2data)

        mergeddata = np.add(wig1data, wig2data)

        merged_wig[chrom] = mergeddata

    return merged_wig


def write_to_wig(wigdata, wigfile):

    with open(wigfile, 'w') as wigout:
        for chrom_key, data in wigdata.items():
            header = "fixedStep chrom={} start={} step={} span={}\n".format(*chrom_key)
            wigout.write(header)
            for binval in data:

                binval = str(binval)

                wigout.write(binval + "\n")


def main(wigs_to_merge, output_wig):
    mergedwigs = None

    for wigfile in wigs_to_merge:
        wigdata = read_wig(wigfile, counts=True)

        if mergedwigs:
            mergedwigs = merge_two_wigs(mergedwigs, wigdata)
        else:
            mergedwigs = wigdata

    write_to_wig(mergedwigs, output_wig)


def parse_args():
    '''
    specify and parse args
    '''

    parser = argparse.ArgumentParser(description='''merge tsv/csv files''')

    parser.add_argument('--input',
                        required=True,
                        nargs="*",
                        help='''wig files to merge''')

    parser.add_argument('--output',
                        required=True,
                        help='''path to merged wig file  ''')

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = parse_args()

    main(args.input, args.output)
