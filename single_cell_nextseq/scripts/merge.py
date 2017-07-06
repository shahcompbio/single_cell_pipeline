'''
Created on Oct 17, 2016

@author: dgrewal
'''
import argparse
from merge_tsv import MergeFiles


def parse_args():
    '''
    specify and parse args
    '''

    parser = argparse.ArgumentParser(description='''merge tsv/csv files''')

    parser.add_argument('--input',
                            nargs='*',
                            required=True,
                            help='''input files for concatenation ''')

    parser.add_argument('--separator',
                            required=True,
                            default="comma",
                            choices = ("comma","tab"),
                            help='''separator type, comma for csv, tab for tsv''')

    parser.add_argument('--type',
                            required=True,
                            default="concatenate",
                            choices = ("concatenate","merge"),
                            help='''merge type''')

    parser.add_argument('--key_cols',
                            nargs='*',
                            help='''keys for merge (column names)''')

    parser.add_argument('--merge_type',
                            default='outer',
                            choices=('inner', 'outer'),
                            help='''method for selection output col names,\
                            outer will keep all columns, inner will keep\
                            columns that exist in both files''')

    parser.add_argument('--nan_value',
                            default='NA',
                            help='value to fill in the missing values')

    parser.add_argument('--output',
                        required=True,
                        help='''path to output file''')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    ARGS = parse_args()
    m = MergeFiles(ARGS)
    m.main()
