# -*- coding: utf-8 -*-
"""
@author: dgrewal

Last updated: Diljot Grewal <dgrewal@bccrc.ca> Jun 3 2015

reads vcf files, filters and write the output in tsv format
"""

#!/usr/bin/env python

import argparse
from vizutils import Utils as pau
from vizutils import Vcf

class ParseMuseq(object):
    '''
    parse, filter and print museq vcf in tsv format
    '''

    def __init__(self, args):
        self.args = args

        self.infiles = pau.get_inputs(self.args.tumour_id,
                                      self.args.normal_id,
                                      self.args.case_id,
                                      self.args.infile,
                                      self.args.all_files,
                                      fh_names='infile')

        pau.test_args(self.args)

        self.genes = pau.read_file_to_list(self.args.genes)

    def main(self):
        '''
        loop through files, load, filter and print
        '''
        header = False
        with open(self.args.output, 'w') as outfile:
            for (case, tum, norm), fname in self.infiles.iteritems():
    
                museq = Vcf(tumour_id = tum,
                            normal_id = norm,
                            case_id = case,
                            project = self.args.project,
                            infile = fname,
                            snpeff_keywords = self.args.snpeff_keywords,
                            keep_dbsnp = self.args.keep_dbsnp,
                            keep_1000gen = self.args.keep_1000gen,
                            chromosomes = self.args.chromosomes,
                            genes = self.genes,
                            rmdups = self.args.remove_duplicates,
                            pr_threshold = self.args.pr_threshold,
                            mode = 'museq'
                           )
                #write header
                if not header:
                    colnames = museq.get_info_header()
                    pau.write_list(outfile, colnames)
                    header=True
                
                infos = museq.get_data()
    
                for info in infos:
                    pau.write_list(outfile, info)


def parse_args():
    '''
    specify and load command line params
    '''
    parser = argparse.ArgumentParser(description='''reads vcf files, filters and\
                                    writes the output in tsv format''')

    exgroup = parser.add_mutually_exclusive_group(required=True)
    exgroup.add_argument("--all_files",
                         default=None,
                         help="""path to the file that has a list of\
                             files and their corresponding ids""")

    exgroup.add_argument("--infile",
                         default=None,
                         help="please provide path to mutationseq's vcf file")

    parser.add_argument("--tumour_id",
                        help='''tumour id for the infile
                            (only required when infile is specified)''')

    parser.add_argument("--normal_id",
                        help='''normal id for the infile
                            (only required when infile is specified)''')

    parser.add_argument("--case_id",
                        help='''case id for the infile
                            (only required when infile is specified)''')

    parser.add_argument("--output",
                        help="Output file name/path")

    parser.add_argument("--snpeff_keywords",
                        default=None,
                        nargs='*',
                        help='''remove calls where the annotation falls in\
                                 this list, default: no filtering''')

    parser.add_argument('--pr_threshold',
                        default=0.5,
                        type=float,
                        help=''' filters calls with probability\
                        below the threshold
                        ''')

    parser.add_argument('--genes',
                        help='''filters all genes except the\
                        ones specified here (default: no filtering)
                        ''')

    parser.add_argument("--project",
                        default='project',
                        help="The project name for the input files")

    parser.add_argument("--keep_dbsnp",
                        action='store_true',
                        default=False,
                        help="The parser won't filter out the "
                        "dbsnp=T calls if flag is set")

    parser.add_argument("--keep_1000gen",
                        action='store_true',
                        default=False,
                        help="The parser won't filter out the "
                        "1000gen=T calls if flag is set")

    parser.add_argument("--keep_cosmic",
                        action='store_true',
                        default=False,
                        help="The parser won't filter out the "
                        "Cosmic=(...) calls if flag is set")

    parser.add_argument("--remove_duplicates",
                        default=False,
                        action='store_true',
                        help=''' Each call will result in multiple lines in\
                                output, one for each gene, default: only the\
                                first highest impact gene will be printed ''')

    parser.add_argument("--chromosomes",
                        default=None,
                        nargs='*',
                        help='''all calls that don't fall in the specified\
                                chromosomes will be filtered''')

    parser.add_argument("--mappability_ref",
                        default=None,
                        help="The reference file with the coordinates of low mappability regions"
                        " These regions will be filtered if the file is provided.")

    parser.add_argument("--rscript_path",
                        default='Rscript',
                        help="path to the Rscript executable")

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    ARGS = parse_args()

    PARSEMUSEQ = ParseMuseq(ARGS)
    PARSEMUSEQ.main()