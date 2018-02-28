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


def _run_museq_worker(tumour, normal, out, log, config, interval):
    '''
    Run museq script for each chromosome

    :param tumour: path to tumour bam
    :param normal: path to normal bam
    :param out: path to temporary output VCF file for the chromosome
    :param log: path to the log file
    :param config: path to the config YAML file
    :param chrom: chromosome number
    '''

    interval = interval.split('_')
    interval = interval[0]+':'+interval[1]+"-"+interval[2]

    script = os.path.join(config['mutationseq'], 'classify.py')
    conf = os.path.join(config['mutationseq'], 'metadata.config')
    model = config['mutationseq_model']
    reference = config['ref_genome']

    cmd = ['python', script, 'normal:' + normal, 'tumour:' + tumour,
           'reference:' + reference, 'model:' + model, '--out', out,
           '--log', log, '--config', conf, '--interval', interval]


    subprocess.call(cmd)


def run_museq(tumour, tumour_bai, normal, normal_bai, out, log, tempdir, config, intervals, ncores=None):
    '''
    Run museq script for all chromosomes and merge VCF files

    :param tumour: path to tumour bam
    :param normal: path to normal bam
    :param out: path to the temporary output VCF file for the merged VCF files
    :param log: path to the log file
    :param config: path to the config YAML file
    '''

    helpers.makedirs(os.path.dirname(log))
    helpers.makedirs(tempdir)

    count = multiprocessing.cpu_count()
    
    if ncores:
        count = min(ncores, count)
    
    pool = multiprocessing.Pool(processes=count)

    output_vcf = {}

    tasks = []

    for interval in intervals:
        outfile = os.path.join(tempdir, str(interval) + '.vcf.tmp')

        task = pool.apply_async(_run_museq_worker,
                                args=(
                                        tumour,
                                        normal,
                                        outfile,
                                        log,
                                        config,
                                        interval,
                                    )
                                )

        tasks.append(task)
        output_vcf[interval] = outfile

    pool.close()
    pool.join()

    [task.get() for task in tasks]

    vcfutils.concatenate_vcf(output_vcf, out)

def parse_museq(infile, output):
    parser = ParseMuseq(infile=infile, tid='NA', nid='NA', output=output,
                        keep_dbsnp=True,keep_1000gen=True,
                        remove_duplicates=True)
    
    parser.main()
