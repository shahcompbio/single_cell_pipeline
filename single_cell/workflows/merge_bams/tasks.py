'''
Created on Jul 24, 2017

@author: dgrewal
'''
import subprocess
import multiprocessing
from single_cell.utils import bamutils
from single_cell.utils import helpers

from subprocess import Popen, PIPE


def merge_bam_worker(input_bam_files, output_bam, output_bai, region):

    #pypeliner.commandline.execute doesnt work with multiprocess.
    #not sure why
#     bamutils.bam_merge(input_bam_files, output_bam, region=region)
#     bamutils.bam_index(output_bam, output_bai)

    cmd = ['samtools', 'merge', '-f', '-R', region]
    cmd.append(output_bam)
    cmd.extend(input_bam_files)
    helpers.run_cmd(cmd)

    cmd = ['samtools',  'index', output_bam, output_bai]
    helpers.run_cmd(cmd)

def merge_bams(bams, bais, outputs, output_index, regions, ncores=None):

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
        task = pool.apply_async(merge_bam_worker,
                         args=(bams, output_bam, output_bai, region)
                        )
        tasks.append(task)

    pool.close()
    pool.join()
    
    [task.get() for task in tasks]

