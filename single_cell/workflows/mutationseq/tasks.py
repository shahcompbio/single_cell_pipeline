'''
Created on Jul 24, 2017

@author: dgrewal
'''
import pypeliner
from single_cell.utils import vcfutils


def subsample(input_bam, output_bam, max_coverage=1000):
    cmd = ['variant', input_bam, '-m', max_coverage, '-v', '-b', '-o', output_bam]
    pypeliner.commandline.execute(*cmd)

    cmd = ['samtools', 'index', output_bam]
    pypeliner.commandline.execute(*cmd)


def run_museq(tumour, normal, out, log, region, config):
    '''
    Run museq script for each chromosome

    :param tumour: path to tumour bam
    :param normal: path to normal bam
    :param out: path to temporary output VCF file for the chromosome
    :param log: path to the log file
    :param config: path to the config YAML file
    :param chrom: chromosome number
    '''

    reference = config['ref_genome']

    region = '{}:{}-{}'.format(*region.split('-'))

    cmd = ['museq', 'normal:' + normal, 'tumour:' + tumour,
           'reference:' + reference, '--out', out,
           '--log', log, '--interval', region]

    museq_params = config.get('museq_params', {})
    for key, val in museq_params.items():
        if isinstance(val, bool):
            if val:
                cmd.append('--{}'.format(key))
        else:
            cmd.append('--{}'.format(key))
            if isinstance(val, list):
                cmd.extend(val)
            else:
                cmd.append(val)

    pypeliner.commandline.execute(*cmd)


def concatenate_vcfs(inputs, output):
    vcfutils.concatenate_vcf(inputs, output)
