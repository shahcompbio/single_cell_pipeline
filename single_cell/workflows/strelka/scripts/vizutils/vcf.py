'''
Created on Apr 10, 2017

@author: dgrewal
'''
from utils import Utils
from parseutils import ParseUtils
from ast import literal_eval as safe_eval
import logging
import gzip
class Vcf(ParseUtils):
    def __init__(self, **kwargs):
        self.tumour_id = kwargs.get("tumour_id")
        self.normal_id = kwargs.get("normal_id")
        self.case_id = kwargs.get("case_id")
        self.project = kwargs.get("project")
        
        self.infile = kwargs.get("infile")
        self.labs = kwargs.get("labels")
        self.labs = self.labs if self.labs else []
        self.lab_vals = kwargs.get("label_vals")
        self.lab_vals = self.lab_vals if self.lab_vals else []


        self.snpeff_keywords = kwargs.get("snpeff_keywords")
        self.keep_dbsnp = kwargs.get("keep_dbsnp")
        self.keep_1000gen = kwargs.get("keep_1000gen")
        self.chromosomes = kwargs.get("chromosomes")
        self.genes = kwargs.get("genes")
        self.rmdups = kwargs.get("rmdups")
        
        self.mode = kwargs.get("mode")
        self.pr_threshold = kwargs.get("pr_threshold")


    @staticmethod
    def is_gzip(filename):
        """
        Uses the file contents to check if the file is gzip or not.
        The magic number for gzip is 1f 8b
        See KRONOS-8 for details
        """
        with open(filename) as f:
            file_start = f.read(4)
        
            if file_start.startswith("\x1f\x8b\x08"):
                return True
            return False
    

    def get_info_header(self):
        '''
        returns header names
        '''

        cols = ['case_id', 'normal_id', 'tumour_id', 'chromosome', 'start',
                'stop', 'gene', 'gene_id', 'type', 'ref', 'alt', 'tr', 'ta',
                'nr', 'na', 'dbsnp', 'thousand_genomes', 'cosmic', 'caller',
                'amino_acid_change', 'functional_class', 'gene_coding',
                'project', 'dna_change', 'biotype', 'filter', 'impact',
                'substitution', 'gt', 'pl', 'mut_pr', 'trinucleotide_ref',
                 'trinucleotide_alt'] + self.labs
        return cols

    def __parse_header_strelka(self, fname):
        '''
        get strelka's mode from header
        i.e. snv or indel
        '''

        if self.is_gzip(fname):
            freader = gzip.open(fname, 'rb')
        else:
            freader = open(fname)

        caller = None
        for line in freader:
            if line[0] != '#':
                break
            line = line.strip().split('=')
            if line[0] == '##content':
                if 'snv' in line[1]:
                    caller = 'strelka_snv'
                elif 'indel' in line[1]:
                    caller = 'strelka_indel'

        if not caller:
            raise Exception('couldn\'t get the information from header')

        return caller

    def parse_strelka(self, fname):
        '''
        extracts all the information about each mutations from line
        '''
        mode = self.__parse_header_strelka(fname)

        if self.is_gzip(fname):
            freader = gzip.open(fname, 'rb')
        else:
            freader = open(fname)

        for line in freader:
            if line[0] == '#':
                continue

            line = line.strip().split('\t')

            line[7] = line[7].split(';')
            line[8] = Utils.build_indices(line[8].split(':'))
            #normal
            line[9] = [safe_eval(val) for val in line[9].split(':')]
            #tumour
            line[10] = [safe_eval(val) for val in line[10].split(':')]

            if mode == 'strelka_indel':
                #TR = (DP+DP2)-sum(TIR)
                tum_alt = line[10][line[8]['TIR']][0]
                tum_ref = line[10][line[8]['DP']] - tum_alt
                norm_alt = line[9][line[8]['TIR']][0]
                norm_ref = line[9][line[8]['DP']] - norm_alt
            else:
                tum_ref = line[10][line[8][line[3]+'U']][0]
                norm_ref = line[9][line[8][line[3]+'U']][0]

                if line[4] == '.':
                    tum_alt = 'N/A'
                    norm_alt = 'N/A'
                else:
                    tum_alt = max([line[10][line[8][val+'U']][0] for val in
                                   line[4].split(',')])

                    norm_alt = max([line[9][line[8][val+'U']][0] for val in
                                    line[4].split(',')])


            subpat, _ = Utils.get_sub_pattern(line[3], line[4])
            #get snpeff annotations
            snpeff = self.parse_snpeff(line[7])
            snpeff = self.sort_snpeff(snpeff)
            #only keep the highest ranked annotation when removing duplicates
            if self.rmdups:
                snpeff = [snpeff[0]]

            anns = self.get_annotations(line[7])

            for snpval in snpeff:
                outval = [self.case_id, self.normal_id, self.tumour_id, line[0],
                          line[1], line[1], snpval[1], snpval[2], snpval[0],
                          line[3], line[4], tum_ref, tum_alt, norm_ref, norm_alt,
                          anns[0], anns[1], anns[2], mode, snpval[4], snpval[5],
                          snpval[6], self.project, snpval[7], snpval[8], line[6],
                          snpval[3], subpat] + [None]*5 +  self.lab_vals


                outval = [val if val else 'N/A' for val in outval]
                yield outval

        freader.close()


    def __parse_header_museq(self, fname):
        '''
        get the bam file paths
        check if museq is single mode
        '''
        header = []

        if self.is_gzip(fname):
            freader = gzip.open(fname, 'rb')
        else:
            freader = open(fname)

        for line in freader:
            if line[0] == '#':
                header.append(line)
                continue
            break
        freader.close()

        for line in header:
            line = line.strip().split('=')
            if line[0] == '##normal':
                nfile = line[1]
            elif line[0] == '##tumour':
                tfile = line[1]
            elif line[0] == '##model':
                if 'single' in line[1] and (tfile == 'N/A' or nfile == 'N/A'):
                    single = True
                else:
                    single = False

        return single

    # pylint: disable=too-many-locals
    # could use a single list but names make code readable
    def parse_museq(self, fname):
        '''
        extracts all the information about each mutations from line
        '''
        single = self.__parse_header_museq(fname)

        caller = 'mutationseq_ss' if single else 'mutationseq'

        if self.is_gzip(fname):
            freader = gzip.open(fname, 'rb')
        else:
            freader = open(fname)

        for line in freader:

            if line[0] == '#':
                continue

            line = line.strip().split()

            chrom = line[0]
            pos = line[1]
            ref = line[3]
            alt = line[4]

            fltr = line[6]

            info = line[7]
            info = info.split(';')
            prob = info[0].split('=')[1]
            tum_ref = info[1].split('=')[1]
            tum_alt = info[2].split('=')[1]
            norm_ref = info[3].split('=')[1]
            norm_alt = info[4].split('=')[1]

            trinuc = info[5].split('=')[1]
            tc_ref = trinuc
            tc_alt = ''.join([trinuc[0], alt, trinuc[2]])

            geno = 'N/A'
            lkl = 'N/A'
            if single:
                try:
                    geno = info[8].split('=')[1]
                    lkl = info[9].split('=')[1]
                except IndexError:
                    logging.getLogger("single_cell.wgs.parse").warn(
                        'Single mode mutationseq but the genotype information is missing'
                    )

            dbsnp, th_gen, cosmic = self.get_annotations(info)

            # get snpeff annotations
            snpeff = self.parse_snpeff(line[7])
            snpeff = self.sort_snpeff(snpeff)

            # only keep the highest ranked annotation when removing duplicates
            if self.rmdups:
                snpeff = [snpeff[0]]

            subpat, _ = Utils.get_sub_pattern(ref, alt, tc_ref)

            for snpval in snpeff:
                (keyword, gene_name, gene_id, impact, amino,
                 func_cl, gene_coding, dna_change, biotype) = snpval
                
                info = [self.case_id, self.normal_id, self.tumour_id,
                        chrom, pos, pos,gene_name, gene_id, keyword, ref,
                        alt, tum_ref, tum_alt, norm_ref, norm_alt, dbsnp, th_gen,
                        cosmic, caller, amino, func_cl, gene_coding, self.project,
                        dna_change, biotype, fltr, impact, subpat, geno, lkl,
                        prob, tc_ref, tc_alt] + self.lab_vals

                info = [val if val else 'N/A' for val in info]
                yield info
        freader.close()

    def filter(self, infos):
        '''
        filter unnecessary calls
        '''
        for info in infos:
            if self.genes and info[6] not in self.genes:
                continue

            if self.snpeff_keywords and \
                info[8] not in self.snpeff_keywords:
                continue

            if info[15] not in ['F', 'N/A'] and not self.keep_dbsnp:
                continue

            if info[16] not in ['F', 'N/A'] and not self.keep_1000gen:
                continue

            if self.chromosomes and info[3] not in self.chromosomes:
                continue

            if self.mode == 'museq' and float(info[30]) < self.pr_threshold:
                continue


            yield info

    def get_data(self):
        """
        get the generator object with the outdata
        """
        if self.mode == 'museq':
            infos = self.parse_museq(self.infile)
        else:
            infos = self.parse_strelka(self.infile)

        infos = self.filter(infos)

        return infos
