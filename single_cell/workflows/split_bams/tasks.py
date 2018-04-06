'''
Created on Nov 21, 2017

@author: dgrewal
'''
import pypeliner

from single_cell.utils import helpers

def split_bam_worker(bam, output_bam, region):

    region = '{}:{}-{}'.format(*region.split('-'))
    cmd = ['samtools', 'view', '-b', bam, '-o', output_bam, region]

    helpers.run_cmd(cmd)

def index_bam_worker(bam, bai):
    ##TODO: specifying bai file path corrupts bams
    cmd = ['samtools', 'index', bam, bai]
    helpers.run_cmd(cmd)


def split_bam_file_one_job(bam, bai, outbam, outbai, regions, ncores=None):

    args = [(bam, outbam(region), region) for region in regions]

    helpers.run_in_parallel(split_bam_worker, args, ncores=ncores)


    args = [(outbam(region), outbai(region)) for region in regions]

    helpers.run_in_parallel(index_bam_worker, args, ncores=ncores)


def split_bam_file(bam, bai, outbam, outbai, interval):

    pypeliner.commandline.execute(
        'samtools', 'view', '-b', bam, interval,
        '>', outbam)

    pypeliner.commandline.execute(
        'samtools', 'index', outbam,
        outbai)
