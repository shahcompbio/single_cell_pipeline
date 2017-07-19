# -*- coding: utf-8 -*-
"""
@author: dgrewal

Last updated: Diljot Grewal <dgrewal@bccrc.ca> Feb 3 2015

This script reads and filters vcf files and generates an output file in tsv format
"""

#!/usr/bin/env python

import argparse
from vizutils import Utils as pau
from vizutils import Vcf


class ParseStrelka(object):
    '''
    parses strelka files
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

        self.index = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

        self.genes = pau.read_file_to_list(self.args.genes)


    def write(self, infos, tum_samp, ref_samp, case, outfile):
        '''
        write info to output.
        Function loops over all data, groups the calls that have same
        pos, gets one with longest alt(longest indel) and prints to file.
        '''
        gr_info = []
        for info in infos:
            if not gr_info:
                gr_info = [info]

            if gr_info[0][3:5] != info[3:5]:

                # get the info line with longest alt length
                lens = [len(val[10]) for val in gr_info]

                out_info = gr_info[lens.index(max(lens))]

                pau.write_list(outfile, out_info)
                #refresh buffer
                gr_info = []

            #add the vals with same pos to list
            gr_info.append(info)

        #write the last value
        if gr_info:
            pau.write_list(outfile, gr_info[0])

    def main(self):
        '''
        loop through files, load, filter and print
        '''

        header = False
        with open(self.args.output, 'w') as outfile:
            for (case, tum_samp, ref_samp), fname in self.infiles.iteritems():
                strelka = Vcf(tumour_id = tum_samp,
                                  normal_id = ref_samp,
                                  case_id = case,
                                  project = self.args.project,
                                  infile = fname,
                                  snpeff_keywords = self.args.snpeff_keywords,
                                  keep_dbsnp = self.args.keep_dbsnp,
                                  keep_1000gen = self.args.keep_1000gen,
                                  chromosomes = self.args.chromosomes,
                                  genes = self.genes,
                                  rmdups = self.args.remove_duplicates
                                  )
    
                #write header
                if not header:
                    colnames = strelka.get_info_header()
                    pau.write_list(outfile, colnames)
                    header=True

    
                infos = strelka.get_data()
    
                self.write(infos, tum_samp, ref_samp, case, outfile)

def parse_args():
    '''
    specify and parse command line args
    '''
    parser = argparse.ArgumentParser(description='''This script reads and
                                    filters vcf files and generates an output
                                    file in tsv format''')

    exgroup = parser.add_mutually_exclusive_group(required=True)
    exgroup.add_argument("--all_files",
                         default=None,
                         help="""path to the file with a list of input files\
                         and their corresponding ids """)

    exgroup.add_argument("--infile",
                         default=None,
                         help="path to strelka's output vcf")

    parser.add_argument("--tumour_id",
                        help='''tumour id for the infile\
                            (only required when infile is specified)''')

    parser.add_argument("--normal_id",
                        help='''normal id for the infile\
                            (only required when infile is specified)''')

    parser.add_argument("--case_id",
                        help='''case id for the infile\
                            (only required when infile is specified)''')

    parser.add_argument("--output",
                        help="output file name")

    parser.add_argument("--genes",
                        default=None,
                        help="remove all genes except the specified")

    parser.add_argument("--snpeff_keywords",
                        default=None,
                        nargs='*',
                        help='''remove calls where the annotation falls in \
                            this list, default: no filtering''')

    parser.add_argument("--project",
                        default='project',
                        help="The project name for the input files")

    parser.add_argument("--remove_duplicates",
                        default=False,
                        action='store_true',
                        help=''' Each call will result in multiple lines in\
                                output, one for each gene, default: only the\
                                first gene will be printed '''
                       )

    parser.add_argument("--keep_dbsnp",
                        action='store_true',
                        default=False,
                        help="The parser won't filter out the dbsnp=T "
                        "calls if flag is set")

    parser.add_argument("--keep_1000gen",
                        action='store_true',
                        default=False,
                        help="The parser won't filter out the 1000gen=T "
                        "calls if flag is set")

    parser.add_argument("--chromosomes",
                        default=None,
                        nargs='*',
                        help='''remove calls that don't fall in the provided chromosomes''')

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

    STRELKA_PARSER = ParseStrelka(ARGS)
    STRELKA_PARSER.main()