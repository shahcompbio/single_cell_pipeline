"""
Created on Feb 19, 2018

@author: dgrewal
"""
import argparse
import json
import os

import pypeliner
from single_cell import __version__
from single_cell.config import generate_batch_config
from single_cell.config import generate_pipeline_config


class parseSentinelPattern(argparse.Action):

    def __call__(self, parser, args, values, option_string=None):
        if values.startswith("/"):
            values = values[1:]

        values = values.split("/")

        parsed = [values[0], "/".join(values[1:])]

        setattr(args, self.dest, parsed)


class parseRegionTemplate(argparse.Action):

    def __call__(self, parser, args, values, option_string=None):
        if "{region}" not in values and "{reads}" not in values:
            raise parser.error(
                "the split template should contain {region} or {reads}"
            )
        setattr(args, self.dest, values)


def separate_pipelinedir_by_subcommand(args):
    if args['which'] in ['generate_config', 'clean_sentinels']:
        return args

    pipelinedir = args.get("pipelinedir", None)

    if pipelinedir:
        subcommand_name = args.get("which", None)

        if not subcommand_name:
            return args

        args["pipelinedir"] = os.path.join(pipelinedir, subcommand_name)

    return args


def add_global_args(parser, dont_add_input_yaml=False):
    pypeliner.app.add_arguments(parser)

    if not dont_add_input_yaml:
        parser.add_argument("--input_yaml",
                            required=True,
                            help='''yaml file with fastq files, output bams and cell metadata''')

    parser.add_argument("--out_dir",
                        required=True,
                        help='''Path to output directory.''')

    config_args = parser.add_mutually_exclusive_group()

    config_args.add_argument(
        "--config_file",
        help='''Path to output directory.''')

    config_args.add_argument(
        "--config_override",
        type=json.loads,
        help='''json string to override the defaults in config''')

    parser.add_argument('--run_with_docker', help='launches pipeline in a docker container')

    return parser


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--version', action='version',
                        version='{version}'.format(version=__version__))

    subparsers = parser.add_subparsers()

    # ===========
    # align
    # ===========
    alignment = add_global_args(subparsers.add_parser("alignment"))
    alignment.set_defaults(which='alignment')
    alignment.add_argument("--library_id",
                           required=True,
                           help='''Library id.''')

    # ===========
    # hmmcopy
    # ===========
    hmmcopy = add_global_args(subparsers.add_parser("hmmcopy"))
    hmmcopy.set_defaults(which='hmmcopy')
    hmmcopy.add_argument("--library_id",
                         required=True,
                         help='''Library id.''')

    # ===========
    # annotation
    # ===========
    annotation = add_global_args(subparsers.add_parser("annotation"))
    annotation.set_defaults(which='annotation')
    annotation.add_argument("--library_id",
                            required=True,
                            help='''Library id.''')

    # ===========
    # merge bams
    # ===========
    merge_bams = add_global_args(subparsers.add_parser("merge_bams"))
    merge_bams.set_defaults(which='merge_bams')

    # ===========
    # split bam
    # ===========
    split_bam = add_global_args(
        subparsers.add_parser("split_bam"),
        dont_add_input_yaml=True)
    split_bam.set_defaults(which='split_bam')

    split_bam.add_argument("--split_bam_template",
                           action=parseRegionTemplate,
                           required=True,
                           help='''template for saving the bams split by region,
                           use {} as place holder for genomic region''')

    split_bam.add_argument("--wgs_bam",
                           required=True,
                           help='''path to the whole genome bam file''')

    # ================
    # variant calling
    # ================
    variant_calling = add_global_args(subparsers.add_parser("variant_calling"))
    variant_calling.set_defaults(which='variant_calling')

    # ===========
    # titan, remixt
    # ===========
    copy_number_calling = add_global_args(
        subparsers.add_parser("copy_number_calling"))
    copy_number_calling.set_defaults(which='copy_number_calling')

    copy_number_calling.add_argument("--clone_id",
                                     required=True,
                                     help='''ID to identify the results''')

    # ===========
    # haplotype blocks
    # ===========
    infer_haps = add_global_args(
        subparsers.add_parser("infer_haps"))
    infer_haps.set_defaults(which='infer_haps')

    infer_haps.add_argument("--normal",
                            default=False,
                            action='store_true',
                            help='''code assumes input is tumour, set this flag to to override''')

    # ===========
    # destruct
    # ===========
    breakpoint_calling = add_global_args(
        subparsers.add_parser("breakpoint_calling"))
    breakpoint_calling.set_defaults(which='breakpoint_calling')

    breakpoint_calling.add_argument("--destruct",
                                    action='store_true',
                                    default=False,
                                    help='''run destruct''')
    breakpoint_calling.add_argument("--lumpy",
                                    action='store_true',
                                    default=False,
                                    help='''run lumpy''')

    # ======================================
    # generates pipeline and batch configs
    # ======================================
    generate_config = subparsers.add_parser("generate_config")
    generate_config.set_defaults(which='generate_config')

    generate_config.add_argument("--pipeline_config",
                                 help='''output yaml file''')

    generate_config.add_argument("--batch_config",
                                 help='''output yaml file''')

    generate_config.add_argument("--config_override",
                                 type=json.loads,
                                 help='''json string to override the defaults in config''')

    # ============================
    # remove tasks from sentinels
    # ============================
    clean_sentinels = subparsers.add_parser("clean_sentinels")
    clean_sentinels.set_defaults(which='clean_sentinels')

    clean_sentinels.add_argument("--mode",
                                 required=True,
                                 choices=['list', 'delete'],
                                 help='''list or delete''')

    clean_sentinels.add_argument("--pattern",
                                 action=parseSentinelPattern,
                                 required=True,
                                 help='''pattern to clean''')

    clean_sentinels.add_argument("--pipelinedir",
                                 required=True,
                                 help='''path to the pipeline dir''')

    args = vars(parser.parse_args())

    # add config paths to global args if needed.
    args = generate_pipeline_config.generate_pipeline_config_in_temp(args)

    args = generate_batch_config.generate_submit_config_in_temp(args)

    # separate pipelinedirs of subcommands
    args = separate_pipelinedir_by_subcommand(args)

    return args
