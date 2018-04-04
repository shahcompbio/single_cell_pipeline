'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner
import warnings
import multiprocessing

from scripts import ParseMuseq

from single_cell.utils import helpers
from single_cell.utils import vcfutils

import subprocess 


def run_museq(tumour, tumour_bai, normal, normal_bai, out, log, config, region_idx, all_regions):
    '''
    Run museq script for each chromosome

    :param tumour: path to tumour bam
    :param normal: path to normal bam
    :param out: path to temporary output VCF file for the chromosome
    :param log: path to the log file
    :param config: path to the config YAML file
    :param chrom: chromosome number
    '''

    region = all_regions[region_idx]

    script = os.path.join(config['mutationseq'], 'classify.py')
    conf = os.path.join(config['mutationseq'], 'metadata.config')
    model = config['mutationseq_model']
    reference = config['ref_genome']

    cmd = ['python', script, 'normal:' + normal, 'tumour:' + tumour,
           'reference:' + reference, 'model:' + model, '--out', out,
           '--log', log, '--config', conf, '--interval', region]

    
    pypeliner.commandline.execute(*cmd)


def concatenate_vcfs(inputs, output):
    vcfutils.concatenate_vcf(inputs, output)


def parse_museq(infile, output):
    parser = ParseMuseq(infile=infile, tid='NA', nid='NA', output=output,
                        keep_dbsnp=True,keep_1000gen=True,
                        remove_duplicates=True)
    
    parser.main()
