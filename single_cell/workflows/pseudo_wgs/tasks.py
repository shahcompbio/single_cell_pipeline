'''
Created on Jul 24, 2017

@author: dgrewal
'''
import subprocess
import multiprocessing
from single_cell.utils import bamutils
from single_cell.utils import helpers


def get_bam_regions(bam_file, split_size, chromosomes):
    helpers.get_bam_regions(bam_file, split_size, chromosomes)


def merge_bam_worker(input_bam_files, output_bam, output_bai, region):

    #pypeliner.commandline.execute doesnt work with multiprocess.
    #not sure why
#     bamutils.bam_merge(input_bam_files, output_bam, region=region)
#     bamutils.bam_index(output_bam, output_bai)

    cmd = ['samtools', 'merge', '-R', region]
    cmd.append(output_bam)
    cmd.extend(input_bam_files)
    subprocess.call(cmd)

    cmd = ['samtools',  'index', output_bam, output_bai]
    subprocess.call(cmd)

def merge_bams(inputs, outputs, output_index, regions, ncores=None):

    count = multiprocessing.cpu_count()
    
    if ncores:
        count = min(ncores, count)
    
    pool = multiprocessing.Pool(processes=count)

    tasks = []

    inputs = inputs.values()

    for region_idx,region in regions.iteritems():
        output_bam = outputs[region_idx]
        output_bai = output_index[region_idx]

        task = pool.apply_async(merge_bam_worker,
                         args=(inputs, output_bam, output_bai, region)
                        )
        tasks.append(task)

    pool.close()
    pool.join()
    
    [task.get() for task in tasks]

