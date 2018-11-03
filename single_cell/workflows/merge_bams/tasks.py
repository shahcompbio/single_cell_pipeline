'''
Created on Jul 24, 2017

@author: dgrewal
'''
import subprocess
import multiprocessing
from single_cell.utils import bamutils
from single_cell.utils import helpers

from subprocess import Popen, PIPE


def cell_region_merge_bams(cell_bams, region_bam, region, kwargs):
    cell_bams = cell_bams.values()
    region = '{}:{}-{}'.format(*region.split('-'))

    bamutils.bam_merge(
        cell_bams, region_bam,
        region=region,
        **kwargs)

    bamutils.bam_index(
        region_bam, region_bam+'.bai',
        **kwargs)


def merge_bam_worker(input_bam_files, output_bam, region, kwargs):

    bamutils.bam_merge(
        input_bam_files, output_bam,
        region=region,
        **kwargs)

    bamutils.bam_index(
        output_bam, output_bam+'.bai',
        **kwargs)


def merge_bams(bams, outputs, regions, docker_config, ncores=None):

    count = multiprocessing.cpu_count()

    if ncores:
        count = min(ncores, count)

    pool = multiprocessing.Pool(processes=count)

    tasks = []

    bams = bams.values()

    for region in regions:
        output_bam = outputs[region]

        region = '{}:{}-{}'.format(*region.split('-'))
        task = pool.apply_async(merge_bam_worker,
                         args=(bams, output_bam, region, docker_config)
                        )
        tasks.append(task)

    pool.close()
    pool.join()

    [task.get() for task in tasks]

