#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Created on Aug 25, 2015

@author: dgrewal
'''

from collections import defaultdict
from intervaltree import IntervalTree
import logging
from utils import Utils

class ParseUtils(object):
    '''
    helper functions for the parsers
    in the bigdata pipeline.
    '''

    @staticmethod
    def get_genome_length(build='grch37_hg19'):
        '''
        return the length of genome (inbp)
        for the build provided
        '''
        if build == 'grch37_hg19':

            # based on chr_lengths from hg19 chromosome (build: GRCh37 released
            # February 27, 2009) please see
            # www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/index.shtml

            genome_length = 3095677412
        else:

            raise Exception('Genome length for build %s not available'
                            % build)

        return genome_length

    @staticmethod
    def close_file(ofile):
        '''
        close file object
        '''
        ofile.close()

    @staticmethod
    def parse_pygene(pygene_line):
        '''
        extracts gene symbols from pygenes annotated column
        'ENSG00000223315,RN7SKP279;ENSG00000222774,
        RN7SKP121;ENSG00000221516,AP002806.1;'
        '''

        genes = []
        pygene_line = pygene_line.split(';')

        pygene_line = [val for val in pygene_line if not val == '']

        if not pygene_line:
            # use N/A when no annotation
            return [('N/A', 'N/A')]

        for gene in pygene_line:
            gene_id = gene.split(',')[0]
            gene_name = gene.split(',')[1]
            genes.append((gene_id, gene_name))

        return genes

    @classmethod
    def build_interval_tree(cls, infile, colnames=None):
        """
        builds an interval tree from the segfile
        """

        itree = defaultdict(IntervalTree)

        if not colnames:
            colnames = {'chr': 'chromosome',
                        'start': 'start',
                        'stop': 'stop',
                        'data': None}

        infile_open = open(infile)

        columns = [colnames['chr'], colnames['start'],
                   colnames['stop']]

        if 'data' in colnames and colnames['data']:
            columns.append(colnames['data'])

        header = infile_open.readline()

        idx = Utils.build_indices(header, columns)

        for line in infile_open:
            line = line.strip().split()

            chrom = line[idx[0]]
            start = int(line[idx[1]])
            end = int(line[idx[2]]) + 1
            data = line[idx[3]] if len(columns) == 4 else None

            itree[chrom].addi(start, end, data)
        return itree

    @staticmethod
    def get_annotations(info):
        '''
        extract dbsnp, 1000gen
        and cosmic anns
        '''
        dbsnp = 'N/A'
        th_gen = 'N/A'
        cosmic = 'N/A'

        if isinstance(info, str):
            info = info.split(';')

        for val in info:
            if 'DBSNP' in val:
                dbsnp = val.split('=')[1]
            elif '1000Gen' in val:
                th_gen = val.split('=')[1]
            elif 'Cosmic' in val:
                cosmic = val.split('=')[1]

        return dbsnp, th_gen, cosmic

    @staticmethod
    def get_functional_change(effect):
        """
        use snpeff effect to predict functional change
        """
        fchange = None
        if effect == ' stop-gained':
            fchange = 'NONSENSE'
        elif effect in ['non-synonymous coding', 'missense_variant']:
            fchange = 'MISSENSE'
        elif effect == ['synonymous coding', 'synonymous_variant']:
            fchange = 'SILENT'
        return fchange

    @staticmethod
    def get_gene_coding(biotype):
        """
        use biotypes to guess the gene coding,
        based on : http://uswest.ensembl.org/Help/Faq?id=468

        retained intron was not mentioned in the ensembl list.
        found the following definition:
        retained_intron:= Alternatively spliced transcript that is believed
        to contain intronic sequence relative to other coding transcripts
        adding it to long_noncoding

        added 'pseudogene', 'processed_pseudogene', 'unitary_pseudogene' 
        to pseudogene

        added protein_coding to protein_coding

        IGV_gene -> IG_V_gene
        TRJ_gene -> TR_J_gene
        IGV_pseudogene -> IG_V_pseudogene
        3prime_overlapping_ncrna -> prime3_overlapping_ncrna

        added TR_V_pseudogene to pseudogene
        added TR_V_gene to protein coding

        """
        protein_coding = ['protein_coding', 'nonsense_mediated_decay',
                          'non_stop_decay', 'IG_C_gene', 'IG_D_gene',
                          'IG_gene', 'IG_J_gene', 'IGLV_gene', 'IG_M_gene',
                          'IG_V_gene', 'IG_Z_gene', 'nontranslating_CDS',
                          'polymorphic_pseudogene', 'TR_C_gene', 'TR_D_gene',
                          'TR_J_gene', 'TR_V_gene']

        pseudogene = ['disrupted_domain', 'IG_C_pseudogene', 'IG_J_pseudogene',
                      'IG_pseudogene', 'IG_V_pseudogene',
                      'processed_pseudogene',
                      'transcribed_processed_pseudogene',
                      'transcribed_unitary_pseudogene',
                      'transcribed_unprocessed_pseudogene',
                      'translated_processed_pseudogene', 'TR_J_pseudogene',
                      'unprocessed_pseudogene', 'pseudogene',
                      'processed_pseudogene', 'unitary_pseudogene',
                      'TR_V_pseudogene']

        long_nc = ['prime3_overlapping_ncrna', 'ambiguous_orf', 'antisense',
                   'antisense_RNA', 'lincRNA', 'ncrna_host',
                   'processed_transcript', 'sense_intronic',
                   'sense_overlapping', 'retained_intron']

        short_nc = ['miRNA', 'miRNA_pseudogene', 'misc_RNA',
                    'miscRNA_pseudogene', 'Mt_rRNA', 'Mt_tRNA',
                    'rRNA', 'scRNA', 'snlRNA', 'snoRNA', 'snRNA',
                    'tRNA', 'tRNA_pseudogene']

        if biotype in protein_coding:
            return 'CODING'
        elif biotype in pseudogene:
            return 'NONCODING'
        elif biotype in long_nc:
            return 'NONCODING'
        elif biotype in short_nc:
            return 'NONCODING'
        elif biotype == '':
            return None
        else:
            logging.getLogger("single_cell.wgs.parse").warn(
                'Cant determine gene_coding of %s' % biotype)
            return None

    @classmethod
    def parse_snpeff(cls, info):
        '''
        extract snpeff annotations from info section
        '''
        output = []

        if isinstance(info, str):
            eff = True if 'EFF=' in info else False
            info = info.split(';')
        else:
            eff = True if 'EFF=' in ' '.join(info) else False

        # OCC-216:  older versions of snpeff dont add any annotations in some
        # cases
        if len([val for val in info if val.split('=')[0] == 'EFF']) == 0:
            return [('N/A',) * 9]

        if eff:
            snpeff_ann = [val for val in info if val.split('=')[0] == 'EFF'][0]
            snpeff_ann = snpeff_ann.split('=')[1].split(',')
            for eff in snpeff_ann:
                keyword = eff.split('(')[0]
                eff = eff.split('(')[1].split(')')[0].split('|')
                gene_name = eff[5]
                biotype = eff[6]
                gene_id = eff[8]
                imp = eff[0]
                amino = eff[3]
                dna_change = None
                func_change = eff[1]
                gene_coding = eff[7]

                outval = (keyword, gene_name, gene_id, imp,
                          amino, func_change, gene_coding,
                          dna_change, biotype)
                outval = tuple(val if val else 'N/A' for val in outval)
                output.append(outval)
        else:
            snpeff_ann = [val for val in info if val.split('=')[0] == 'ANN'][0]
            snpeff_ann = snpeff_ann.split('=')[1].split(',')
            for eff in snpeff_ann:
                eff = eff.split('|')
                keyword = eff[1]
                gene_name = eff[3]
                gene_id = eff[4]
                biotype = eff[7]
                imp = eff[2]
                amino = eff[10]
                dna_change = eff[9]
                func_change = cls.get_functional_change(keyword)
                gene_coding = cls.get_gene_coding(biotype)

                outval = (keyword, gene_name, gene_id, imp,
                          amino, func_change, gene_coding,
                          dna_change, biotype)
                outval = tuple(val if val else 'N/A' for val in outval)
                output.append(outval)
        return output

    @staticmethod
    def sort_snpeff(snpeff_ann):
        '''
        sort snpeff annotations by their modifier
        '''
        high = []
        mod = []
        low = []
        oth = []

        for val in snpeff_ann:
            if val[3] == 'HIGH':
                high.append(val)
            elif val[3] == 'MODERATE':
                mod.append(val)
            elif val[3] == 'LOW':
                low.append(val)
            else:
                oth.append(val)

        return high + mod + low + oth

