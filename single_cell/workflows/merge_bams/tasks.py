'''
Created on Jul 24, 2017

@author: dgrewal
'''
import subprocess
import multiprocessing
from single_cell.utils import bamutils
from single_cell.utils import helpers

from subprocess import Popen, PIPE


def merge_bam_worker(input_bam_files, output_bam, output_bai, region, config):

    bamutils.bam_merge(
        input_bam_files, output_bam,
        region=region,
        dockerize=config['dockerize'],
        mounts=config['mounts'],
        image=config['images']['samtools']['image'],
        username=config['images']['samtools']['username'],
        password=config['images']['samtools']['password'],
        server=config['images']['samtools']['server'],
    )

    bamutils.bam_index(
        output_bam, output_bai,
        dockerize=config['dockerize'],
        mounts=config['mounts'],
        image=config['images']['samtools']['image'],
        username=config['images']['samtools']['username'],
        password=config['images']['samtools']['password'],
        server=config['images']['samtools']['server'],
    )


def merge_bams(bams, bais, outputs, output_index, regions, docker_config, ncores=None):

    count = multiprocessing.cpu_count()
    
    if ncores:
        count = min(ncores, count)
    
    pool = multiprocessing.Pool(processes=count)

    tasks = []

    bams = bams.values()

    for region in regions:
        output_bam = outputs[region]
        output_bai = output_index[region]

        region = '{}:{}-{}'.format(*region.split('-'))
        # merge_bam_worker(bams, output_bam, output_bai, region, docker_config)

        task = pool.apply_async(merge_bam_worker,
                         args=(bams, output_bam, output_bai, region, docker_config)
                        )
        tasks.append(task)

    pool.close()
    pool.join()

    [task.get() for task in tasks]

