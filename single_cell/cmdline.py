'''
Created on Feb 19, 2018

@author: dgrewal
'''
import argparse
import pypeliner


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    pypeliner.app.add_arguments(parser)

    subparsers = parser.add_subparsers()

    #Align subcommand
    align = subparsers.add_parser("align")
    align.set_defaults(which='align')


    align.add_argument('sample_info',
                        help='''Per sample meta data CSV''')

    align.add_argument('fastqs_file',
                        help='''Path to input fastq table CSV.''')

    align.add_argument('bams_file',
                        help='''Path to input fastq table CSV.''')

    align.add_argument('library_id',
                        help='''Library id.''')

    align.add_argument('--realign',
                        action='store_true',
                        help='''Lanes to analyze.''')
    

    #hmmcopy command
    hmmcopy = subparsers.add_parser("hmmcopy")
    hmmcopy.set_defaults(which='hmmcopy')

    hmmcopy.add_argument('sample_info',
                        help='''Per sample meta data CSV''')
    hmmcopy.add_argument('bams_file',
                        help='''Path to input fastq table CSV.''')
    hmmcopy.add_argument('library_id',
                        help='''Library id.''')


    #aneufinder command
    aneufinder = subparsers.add_parser("aneufinder")
    aneufinder.set_defaults(which='aneufinder')

    aneufinder.add_argument('sample_info',
                        help='''Per sample meta data CSV''')
    aneufinder.add_argument('bams_file',
                        help='''Path to input fastq table CSV.''')
    aneufinder.add_argument('library_id',
                        help='''Library id.''')


    #hmmcopy command
    summary = subparsers.add_parser("summary")
    summary.set_defaults(which='summary')

    summary.add_argument('sample_info',
                        help='''Per sample meta data CSV''')
    summary.add_argument('library_id',
                        help='''Library id.''')


    #hmmcopy command
    pseudo_wgs = subparsers.add_parser("pseudo_wgs")
    pseudo_wgs.set_defaults(which='pseudo_wgs')

    pseudo_wgs.add_argument('bams_file',
                        help='''Per sample meta data CSV''')


    
    #generate wgs bam command
    #variant calling command
    varcall = subparsers.add_parser("variant_calling")
    varcall.set_defaults(which='variant_calling')

    varcall.add_argument('matched_normal',
                        help='''Path to matched wgs normal.''')
    varcall.add_argument('bams_file',
                        help='''Per sample meta data CSV''')


    # common options
    parser.add_argument('out_dir',
                        help='''Path to output files.''')

    parser.add_argument('config_file',
                        help='''Path to yaml config file.''')

    parser.add_argument('--generate_pseudo_wgs',
                        action='store_true',
                        help='''Lanes to analyze.''')

    args = vars(parser.parse_args())

    return args
