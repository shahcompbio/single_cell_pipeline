'''
Created on Feb 19, 2018

@author: dgrewal
'''
import argparse
import pypeliner


def arg_exists(args, key, mode):
    error_string = "option {} is required in {} mode"

    if key not in args or not args[key]:

        raise Exception(error_string.format(key, mode))
    

def check_update_args(args):

    if args["modes"] == ["qc"]:
        args["modes"] = ["align","hmmcopy","aneufinder","copyclone"]

    if args["modes"] == ["all"]:
        args["modes"] = ["align", "hmmcopy", "aneufinder", "copyclone",\
                         "pseudo_wgs", "split_normal", "variant_calling", "germline_calling"]

    args["merged_wgs_template"] = args["merged_wgs_template"].replace("{}","{region}")
    args["normal_split_template"] = args["normal_split_template"].replace("{}","{region}")

    for mode in args["modes"]:
        if mode in ["align","hmmcopy","aneufinder", "copyclone"]:
            arg_exists(args, "input_yaml", mode)
            arg_exists(args, "library_id", mode)

        elif mode in ["pseudo_wgs"]:
            arg_exists(args, "input_yaml", mode)
            arg_exists(args, "merged_wgs_template", mode)

        elif mode in ["variant_calling"]:
            arg_exists(args, "input_yaml", mode)
            arg_exists(args, "merged_wgs_template", mode)
            arg_exists(args, "normal_split_template", mode)

        elif mode == "split_normal":
            if not args["matched_normal"] and not args["normal_yaml"]:
                raise Exception("matched_normal or normal_yaml required to run split_normal workflow")
            arg_exists(args, "normal_split_template", mode)

        elif mode in ["germline_calling"]:
            arg_exists(args, "input_yaml", mode)
            arg_exists(args, "normal_split_template", mode)
        else:
            raise Exception("unknown mode: {}".format(mode))

    return args

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

    normal_input = parser.add_mutually_exclusive_group()

    normal_input.add_argument('--matched_normal',
                        help='''Path to matched wgs normal.''')

    normal_input.add_argument('--normal_yaml',
                        help='''Path to matched wgs normal.''')

    parser.add_argument('--normal_split_template',
                        help='''Path to matched wgs normal.''')

    parser.add_argument('--merged_wgs_template',
                        help='''Path to pseudo whole genome bam file''')

    parser.add_argument('--realign',
                        action='store_true',
                        help='''will run local realignment on all cells in batch mode''')


    args = vars(parser.parse_args())


    args["modes"] = args["modes"].split(',')


    args = check_update_args(args)

    return args
