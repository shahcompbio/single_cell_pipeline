'''
Created on Feb 19, 2018

@author: dgrewal
'''
import argparse
import pypeliner
import copy
import os
import sys
import json

class parseSentinelPattern(argparse.Action):

    def __call__(self, parser, args, values, option_string=None):
        if values.startswith("/"):
            values = values[1:]

        values = values.split("/")

        parsed = [values[0], "/".join(values[1:])]
    
        setattr(args, self.dest, parsed)


class parseRegionTemplate(argparse.Action):

    def __call__(self, parser, args, values, option_string=None):
        values = values.replace("{}", "{region}")
        setattr(args, self.dest, values)

def arg_exists(args, key, mode, parser):
    error_string = "option {} is required in {} mode"

    if key not in args:
        raise parser.error(error_string.format(key, mode))

def check_required_args(allcommands, parser):

    for mode, args in allcommands.iteritems():

        if mode in ["align","hmmcopy","aneufinder", "copyclone"]:
            arg_exists(args, "input_yaml", mode, parser)
            arg_exists(args, "library_id", mode, parser)

        elif mode in ["merge_bams"]:
            arg_exists(args, "input_yaml", mode, parser)
            arg_exists(args, "merged_bam_template", mode, parser)

        elif mode in ["variant_calling"]:
            arg_exists(args, "input_yaml", mode, parser)
            arg_exists(args, "tumour_template", mode, parser)
            arg_exists(args, "normal_template", mode, parser)

        elif mode in ["split_bam"]:
            arg_exists(args, "split_bam_template", mode, parser)
            arg_exists(args, "wgs_bam", mode, parser)

        elif mode in ["germline_calling"]:
            arg_exists(args, "input_yaml", mode, parser)
            arg_exists(args, "input_template", mode, parser)

        elif mode in ["generate_config", "clean_sentinels"]:
            continue
        else:
            raise parser.error("unknown mode: {}".format(mode))

def generate_args_by_mode(global_args, commands):
    global_args = vars(global_args)

    
    mode_based_args = {}

    for subcommand in commands:
        subcommand = vars(subcommand)
        
        mode = subcommand["which"]
        subcommand.update(global_args)

        if subcommand.get("pipelinedir", None):
            subcommand["pipelinedir"] = os.path.join(subcommand["pipelinedir"], mode) 
        
        mode_based_args[mode] = subcommand

    return mode_based_args


def print_help(globals, subcommands):
    args = sys.argv[1:]

    if "-h" in args or "--help" in args:
        try:
            print subcommands.parse_known_args(sys.argv[1:])
        except:
            print globals.parse_known_args(sys.argv[1:])

def parse_all_commands(globals, subcommands):
    globalargs, unknowns = globals.parse_known_args()

    i=0
    subcommand_args = []
    while unknowns and i<20:
        args, unknowns = subcommands.parse_known_args(unknowns)
        subcommand_args.append(args)
 
    if unknowns:
        parser.error("unrecognized arguments: {}".format(str(unknowns)))
 
    return globalargs, subcommand_args



def parse_args():
    subcommands = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = subcommands.add_subparsers()

    # subparser to align bams
    align = subparsers.add_parser("align")
    align.set_defaults(which='align')
    align.add_argument('--realign',
                       action='store_true',
                       help='''will run local realignment on all cells in batch mode''')

    # subparser to run hmmcopy
    hmmcopy = subparsers.add_parser("hmmcopy")
    hmmcopy.set_defaults(which='hmmcopy')

    # subparser to align bams
    copyclone = subparsers.add_parser("copyclone")
    copyclone.set_defaults(which='copyclone')

    # subparser to align bams
    aneufinder = subparsers.add_parser("aneufinder")
    aneufinder.set_defaults(which='aneufinder')

    # subparser to align bams
    merge_bams = subparsers.add_parser("merge_bams")
    merge_bams.set_defaults(which='merge_bams')
    merge_bams.add_argument("--merged_bam_template",
                            required=True,
                            help='''template for saving the bams merged by region, use {} as place holder for genomic region''')

    # subparser to align bams
    split_bam = subparsers.add_parser("split_bam")
    split_bam.set_defaults(which='split_bam')

    split_bam.add_argument("--split_bam_template",
                           action=parseRegionTemplate,
                           required=True,
                           help='''template for saving the bams merged by region, use {} as place holder for genomic region''')

    split_bam.add_argument("--wgs_bam",
                           required=True,
                           help='''path to the whole genome bam file''')

    # subparser to align bams
    variant_calling = subparsers.add_parser("variant_calling")
    variant_calling.set_defaults(which='variant_calling')

    variant_calling.add_argument("--tumour_template",
                                 required=True,
                                 action=parseRegionTemplate,
                                 help='''template for saving the bams merged by region, use {} as place holder for genomic region''')

    variant_calling.add_argument("--normal_template",
                                 required=True,
                                 action=parseRegionTemplate,
                                 help='''template for saving the bams merged by region, use {} as place holder for genomic region''')

    # subparser to align bams
    germline_calling = subparsers.add_parser("germline_calling")
    germline_calling.set_defaults(which='germline_calling')

    germline_calling.add_argument("--input_template",
                                  required=True,
                                  action=parseRegionTemplate,
                                  help='''template for saving the bams merged by region, use {} as place holder for genomic region''')


    # subparser to align bams
    generate_config = subparsers.add_parser("generate_config")
    generate_config.set_defaults(which='generate_config')

    generate_config.add_argument("--output",
                                  required=True,
                                  help='''output yaml file''')


    clean_sentinels = subparsers.add_parser("clean_sentinels")
    clean_sentinels.set_defaults(which='clean_sentinels')

    clean_sentinels.add_argument("--mode",
                                 required=True,
                                 choices=['list','delete'],
                                 help='''list or delete''')

    clean_sentinels.add_argument("--pattern",
                                 action=parseSentinelPattern,
                                  required=True,
                                  help='''pattern to clean''')


    globals = argparse.ArgumentParser()
    pypeliner.app.add_arguments(globals)

    globals.add_argument("--input_yaml",
                         help='''yaml file with fastq files, output bams and cell metadata''')

    globals.add_argument("--library_id",
                         help='''Library id.''')

    globals.add_argument("--out_dir",
                         help='''Path to output directory.''')

    globals.add_argument("--config_file",
                         help='''Path to output directory.''')

    globals.add_argument("--config_override",
                         type=json.loads,
                         help='''json string to override the defaults in config''')


    #check if -h in args
    print_help(globals, subcommands)

    globalargs, subcommand_args = parse_all_commands(globals, subcommands)

    args = generate_args_by_mode(globalargs, subcommand_args)

    check_required_args(args, subcommands)

    return args

