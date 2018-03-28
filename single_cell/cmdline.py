'''
Created on Feb 19, 2018

@author: dgrewal
'''
import argparse
import pypeliner


def check_args(args):

    for mode in args["modes"]:
        if mode in ["align","hmmcopy","aneufinder", "copyclone"]:
            assert args["input_yaml"]
            assert args["library_id"]

        elif mode in ["pseudo_wgs"]:
            assert args["input_yaml"]
            assert args["merged_wgs"]

        elif mode in ["variant_calling"]:
            assert args["matched_normal"]
            assert args["input_yaml"]
            assert args["merged_wgs"]

        else:
            raise Exception()


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    pypeliner.app.add_arguments(parser)


    parser.add_argument('modes',
                        help='''modes to run, separated by comma.''')

    # common options
    parser.add_argument('out_dir',
                        help='''Path to output files.''')

    parser.add_argument('config_file',
                        help='''Path to yaml config file.''')

    parser.add_argument('--input_yaml',
                        help='''yaml file with fastq files, output bams and cell metadata''')

    parser.add_argument('--library_id',
                        help='''Library id.''')

    parser.add_argument('--matched_normal',
                        help='''Path to matched wgs normal.''')

    parser.add_argument('--merged_wgs',
                        help='''Path to pseudo whole genome bam file''')

    parser.add_argument('--realign',
                        action='store_true',
                        help='''will run local realignment on all cells in batch mode''')


    args = vars(parser.parse_args())


    args["modes"] = args["modes"].split(',')


    check_args(args)

    return args
