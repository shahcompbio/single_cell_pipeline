'''
Created on Feb 19, 2018

@author: dgrewal
'''
import argparse
import pypeliner
import copy
import os

QC_MODES = ["align","hmmcopy","aneufinder","copyclone"]
WGS_MODES = ["pseudo_wgs", "split_normal", "variant_calling", "germline_calling"]

class parseModes(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
        values = values.split(',')
        if "qc" in values:
            values.extend(QC_MODES)
            values.remove("qc")
    
        if "wgs" in values:
            values.extend(WGS_MODES)
            values.remove("wgs")
    
        if "all" in values:
            values.extend(QC_MODES)
            values.extend(WGS_MODES)
            values.remove("all")

        setattr(args, self.dest, values)


class parseRegionTemplate(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
        values = values.replace("{}","{region}")
        setattr(args, self.dest, values)



def arg_exists(args, key, mode, parser):
    error_string = "option {} is required in {} mode"

    if key not in args:
        raise parser.error(error_string.format(key, mode))
    

def check_required_args(parser, args):

    if args["which"] != "run":
        return

    for mode in args["modes"]:

        if mode in ["align","hmmcopy","aneufinder", "copyclone"]:
            arg_exists(args, "input_yaml", mode, parser)
            arg_exists(args, "library_id", mode, parser)

        elif mode in ["pseudo_wgs"]:
            arg_exists(args, "input_yaml", mode, parser)
            arg_exists(args, "merged_wgs_template", mode, parser)

        elif mode in ["variant_calling"]:
            arg_exists(args, "input_yaml", mode, parser)
            arg_exists(args, "merged_wgs_template", mode, parser)
            arg_exists(args, "normal_split_template", mode, parser)

        elif mode in ["split_normal"]:
            if not args["matched_normal"] and not args["normal_yaml"]:
                parser.error("matched_normal or normal_yaml required to run split_normal workflow")
            arg_exists(args, "normal_split_template", mode, parser)

        elif mode in ["germline_calling"]:
            arg_exists(args, "input_yaml", mode, parser)
            arg_exists(args, "normal_split_template", mode, parser)
        else:
            raise parser.error("unknown mode: {}".format(mode))


def generate_args_by_mode(args):

    mode_based_args = {}

    if args["which"] == "generate_config":
        mode_based_args = {"generate_config": args}
    
    elif args["which"] == "run":
        for mode in args["modes"]:
            args_mode = copy.deepcopy(args)
    
            args_mode["pipelinedir"] = os.path.join(args_mode["pipelinedir"], mode)
    
            mode_based_args[mode] = args_mode

    return mode_based_args


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    pypeliner.app.add_arguments(parser)

    subparsers = parser.add_subparsers()

    #All subcommands
    run = subparsers.add_parser("run")
    run.set_defaults(which='run')


    run.add_argument('modes',
                     action=parseModes,
                     help='''modes to run, separated by comma.''')

    # common options
    run.add_argument('out_dir',
                     help='''Path to output files.''')

    run.add_argument('config_file',
                     help='''Path to yaml config file.''')

    run.add_argument('--input_yaml',
                     help='''yaml file with fastq files, output bams and cell metadata''')

    run.add_argument('--library_id',
                     help='''Library id.''')

    normal_input = run.add_mutually_exclusive_group()

    normal_input.add_argument('--matched_normal',
                              help='''Path to matched wgs normal.''')

    normal_input.add_argument('--normal_yaml',
                              help='''Path to matched wgs normal.''')

    run.add_argument('--normal_split_template',
                     action=parseRegionTemplate,
                     help='''Path to matched wgs normal.''')

    run.add_argument('--merged_wgs_template',
                     action=parseRegionTemplate,
                     help='''Path to pseudo whole genome bam file''')

    run.add_argument('--realign',
                     action='store_true',
                     help='''will run local realignment on all cells in batch mode''')


    #All subcommands
    generate_config = subparsers.add_parser("generate_config")
    generate_config.set_defaults(which='generate_config')

    generate_config.add_argument("config_path",
                                 help="path to the config file")


    args = vars(parser.parse_args())

    check_required_args(parser, args)

    args = generate_args_by_mode(args)

    return args
