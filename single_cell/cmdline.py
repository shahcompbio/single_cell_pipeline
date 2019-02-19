"""
Created on Feb 19, 2018

@author: dgrewal
"""
import os
import json
import argparse
import pypeliner
from single_cell.config import generate_batch_config
from single_cell.config import generate_pipeline_config
from single_cell import __version__

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

    parser.add_argument("--config_file",
                        help='''Path to output directory.''')

    parser.add_argument("--config_override",
                        type=json.loads,
                        help='''json string to override the defaults in config''')

    parser.add_argument('--run_with_docker')

    return parser


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--version', action='version',
                        version='{version}'.format(version=__version__))

    subparsers = parser.add_subparsers()

    #===========
    # align
    #===========
    align = add_global_args(subparsers.add_parser("align"))
    align.set_defaults(which='align')

    align.add_argument('--realign',
                       action='store_true',
                       help='''will run local realignment on all cells in batch mode''')

    align.add_argument("--library_id",
                       required=True,
                       help='''Library id.''')

    align.add_argument("--metrics_only",
                       default=False,
                       action="store_true",
                       help='''only run alignment metrics extraction, assumes that bams are available''')


    #===========
    # hmmcopy
    #===========
    hmmcopy = add_global_args(subparsers.add_parser("hmmcopy"))
    hmmcopy.set_defaults(which='hmmcopy')

    hmmcopy.add_argument("--library_id",
                         required=True,
                         help='''Library id.''')

    hmmcopy.add_argument("--alignment_metrics",
                         help='''alignment metrics h5 file''')

    #===========
    # copyclone
    #===========
    copyclone = add_global_args(subparsers.add_parser("copyclone"))
    copyclone.set_defaults(which='copyclone')
    copyclone.add_argument("--library_id",
                           required=True,
                           help='''Library id.''')

    #===========
    # aneufinder
    #===========
    aneufinder = add_global_args(subparsers.add_parser("aneufinder"))
    aneufinder.set_defaults(which='aneufinder')
    aneufinder.add_argument("--library_id",
                            required=True,
                            help='''Library id.''')

    #===========
    # merge bams
    #===========
    merge_bams = add_global_args(subparsers.add_parser("merge_bams"))
    merge_bams.set_defaults(which='merge_bams')
    merge_bams.add_argument("--merged_bam_template",
                            required=True,
                            help='''template for saving the bams merged by region,
                            use {} as place holder for genomic region''')

    #===========
    # split bam
    #===========
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

    #================
    # variant calling
    #================
    variant_calling = add_global_args(subparsers.add_parser("variant_calling"))
    variant_calling.set_defaults(which='variant_calling')

    variant_calling.add_argument("--tumour_template",
                                 required=True,
                                 action=parseRegionTemplate,
                                 help='''template for bams merged by region,
                                 use {} as place holder for genomic region''')

    variant_calling.add_argument("--normal_template",
                                 required=True,
                                 action=parseRegionTemplate,
                                 help='''template for bams merged by region,
                                 use {} as place holder for genomic region''')

    #===========
    # titan, remixt
    #===========
    copy_number_calling = add_global_args(
        subparsers.add_parser("copy_number_calling"), dont_add_input_yaml=True)
    copy_number_calling.set_defaults(which='copy_number_calling')

    copy_number_calling.add_argument("--tumour_yaml",
                                     required=True,
                                     help='''template for bams merged by region,
                                     use {} as place holder for genomic region''')

    copy_number_calling.add_argument("--normal_template",
                                     required=True,
                                     help='''template for bams merged by region,
                                     use {} as place holder for genomic region''')

    copy_number_calling.add_argument("--clone_id",
                                     required=True,
                                     help='''ID to identify the results''')

    #===========
    # haplotype blocks
    #===========
    infer_haps = add_global_args(
        subparsers.add_parser("infer_haps"), dont_add_input_yaml=True)
    infer_haps.set_defaults(which='infer_haps')

    mutexgroup = infer_haps.add_mutually_exclusive_group(required=True)

    mutexgroup.add_argument("--input_yaml",
                            help='''template for bams merged by region,
                                    use {} as place holder for genomic region''')

    mutexgroup.add_argument("--input_bam",
                            help='''template for bams merged by region,
                                    use {} as place holder for genomic region''')

    infer_haps.add_argument("--normal",
                            default=False,
                            action='store_true',
                            help='''code assumes input is tumour, set this flag to to override''')


    #===========
    # germline
    #===========
    germline_calling = add_global_args(
        subparsers.add_parser("germline_calling"))
    germline_calling.set_defaults(which='germline_calling')

    germline_calling.add_argument("--input_template",
                                  required=True,
                                  help='''template for bams merged by region,
                                  use {} as place holder for genomic region''')
    #===========
    # destruct
    #===========
    breakpoint_calling = add_global_args(
        subparsers.add_parser("breakpoint_calling"))
    breakpoint_calling.set_defaults(which='breakpoint_calling')

    breakpoint_calling.add_argument("--matched_normal",
                                    required=True,
                                    help='''normal bam file''')
    breakpoint_calling.add_argument("--normal_id",
                                    required=True,
                                    help='''normal id''')
    breakpoint_calling.add_argument("--tumour_id",
                                    required=True,
                                    help='''id for tumour''')
    breakpoint_calling.add_argument("--destruct",
                                    action='store_true',
                                    default=False,
                                    help='''run destruct''')
    breakpoint_calling.add_argument("--lumpy",
                                    action='store_true',
                                    default=False,
                                    help='''run lumpy''')


    #======================================
    # count variants from multiple samples
    #======================================
    variant_counting = add_global_args(
        subparsers.add_parser("variant_counting"))
    variant_counting.set_defaults(which='variant_counting')

    variant_counting.add_argument("--input_vcfs",
                                  required=True,
                                  help='''vcf files''',
                                  nargs='+')

    #======================================
    # bulk analysis from multiple samples
    #======================================
    multi_sample_pseudo_bulk = add_global_args(
        subparsers.add_parser("multi_sample_pseudo_bulk"))
    multi_sample_pseudo_bulk.set_defaults(which='multi_sample_pseudo_bulk')

    multi_sample_pseudo_bulk.add_argument(
        "--destruct",
        default=False,
        action='store_true',
        help='''add destruct to list of analyses'''
    )


    #======================================
    # generates pipeline and batch configs
    #======================================
    generate_config = subparsers.add_parser("generate_config")
    generate_config.set_defaults(which='generate_config')

    generate_config.add_argument("--pipeline_config",
                                 help='''output yaml file''')

    generate_config.add_argument("--batch_config",
                                 help='''output yaml file''')

    generate_config.add_argument("--config_override",
                                 type=json.loads,
                                 help='''json string to override the defaults in config''')


    #============================
    # remove tasks from sentinels
    #============================
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

    #============================
    # generate LTM tree
    #============================
    ltm = add_global_args(subparsers.add_parser("ltm"), dont_add_input_yaml=True)
    ltm.set_defaults(which='ltm')

    ltm.add_argument("--input_csv",
                     required=True,
                     help='''csv file with alignment and hmmcopy results for each timepoint''')

    ltm.add_argument("--root_id",
                     help='''ID of the cell to use as root of the tree. Default: first SA928 cell''',
                     default=None)

    ltm.add_argument('--number_of_jobs',
                     help='''Number of jobs to submit for distance calculation for scaled method.''',
                     default=10,
                     type=int)

    ltm.add_argument('--ploidy',
                     help='''Ploidy to use for analysis.''',
                     default=0,
                     type=int)


    args = vars(parser.parse_args())

    # add config paths to global args if needed.
    args = generate_pipeline_config.generate_pipeline_config_in_temp(args)

    args = generate_batch_config.generate_submit_config_in_temp(args)

    # separate pipelinedirs of subcommands
    args = separate_pipelinedir_by_subcommand(args)

    return args
