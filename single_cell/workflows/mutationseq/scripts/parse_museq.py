# -*- coding: utf-8 -*-
"""
@author: dgrewal

Last updated: Diljot Grewal <dgrewal@bccrc.ca> Jun 3 2015

reads vcf files, filters and write the output in tsv format
"""

#!/usr/bin/env python

from vizutils import Utils as pau
from vizutils import Vcf

class ParseMuseq(object):
    '''
    parse, filter and print museq vcf in tsv format
    '''

    def __init__(self, **kwargs):

        self.infiles = pau.get_inputs(kwargs.get('tid'),
                                      kwargs.get('nid'),
                                      kwargs.get('case'),
                                      kwargs.get('infile'),
                                      kwargs.get('all_files'),
                                      fh_names='infile')

        self.output = kwargs.get('output')
        self.project = kwargs.get('project')

        self.genes = pau.read_file_to_list(kwargs.get('genes'))
        self.snpeff_keywords = kwargs.get('snpeff_keywords')
        self.chromosomes = kwargs.get('chromosomes')
        self.remove_duplicates = kwargs.get('rm_dups')
        self.pr_threshold = kwargs.get('pr_thres')

        self.keep_dbsnp = kwargs.get('keep_dbsnp')
        self.keep_1000gen = kwargs.get('keep_1000gen')


    def main(self):
        '''
        loop through files, load, filter and print
        '''
        header = False
        with open(self.output, 'w') as outfile:
            for (case, tum, norm), fname in self.infiles.items():
    
                museq = Vcf(tumour_id = tum,
                            normal_id = norm,
                            case_id = case,
                            infile = fname,
                            snpeff_keywords = self.snpeff_keywords,
                            keep_dbsnp = self.keep_dbsnp,
                            keep_1000gen = self.keep_1000gen,
                            chromosomes = self.chromosomes,
                            genes = self.genes,
                            rmdups = self.remove_duplicates,
                            pr_threshold = self.pr_threshold,
                            mode = 'museq'
                           )
                #write header
                if not header:
                    colnames = museq.get_info_header()
                    pau.write_list(outfile, colnames, sep=",")
                    header=True
                
                infos = museq.get_data()
    
                for info in infos:
                    pau.write_list(outfile, info, sep=',')
