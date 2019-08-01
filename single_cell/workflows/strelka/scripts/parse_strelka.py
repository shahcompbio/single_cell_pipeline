# -*- coding: utf-8 -*-
"""
@author: dgrewal

Last updated: Diljot Grewal <dgrewal@bccrc.ca> Feb 3 2015

This script reads and filters vcf files and generates an output file in tsv format
"""

#!/usr/bin/env python

from vizutils import Utils as pau
from vizutils import Vcf

class ParseStrelka(object):
    '''
    parses strelka files
    '''
    def __init__(self, **kwargs):

        self.infiles = pau.get_inputs(kwargs.get('tid'),
                                      kwargs.get('nid'),
                                      kwargs.get('case'),
                                      kwargs.get('infile'),
                                      kwargs.get('all_files'),
                                      fh_names='infile')

        self.index = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

        self.output = kwargs.get('output')
        self.project = kwargs.get('project')

        self.genes = pau.read_file_to_list(kwargs.get('genes'))
        self.snpeff_keywords = kwargs.get('snpeff_keywords')
        self.chromosomes = kwargs.get('chromosomes')
        self.remove_duplicates = kwargs.get('rm_dups')
        self.pr_threshold = kwargs.get('pr_thres')

        self.keep_dbsnp = kwargs.get('keep_dbsnp')
        self.keep_1000gen = kwargs.get('keep_1000gen')

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

                pau.write_list(outfile, out_info, sep=",")
                #refresh buffer
                gr_info = []

            #add the vals with same pos to list
            gr_info.append(info)

        #write the last value
        if gr_info:
            pau.write_list(outfile, gr_info[0], sep=",")

    def main(self):
        '''
        loop through files, load, filter and print
        '''

        header = False
        with open(self.output, 'w') as outfile:
            for (case, tum_samp, ref_samp), fname in self.infiles.items():
                strelka = Vcf(tumour_id = tum_samp,
                                  normal_id = ref_samp,
                                  case_id = case,
                                  project = self.project,
                                  infile = fname,
                                  snpeff_keywords = self.snpeff_keywords,
                                  keep_dbsnp = self.keep_dbsnp,
                                  keep_1000gen = self.keep_1000gen,
                                  chromosomes = self.chromosomes,
                                  genes = self.genes,
                                  rmdups = self.remove_duplicates
                                  )
    
                #write header
                if not header:
                    colnames = strelka.get_info_header()
                    pau.write_list(outfile, colnames, sep=',')
                    header=True

    
                infos = strelka.get_data()
    
                self.write(infos, tum_samp, ref_samp, case, outfile)

