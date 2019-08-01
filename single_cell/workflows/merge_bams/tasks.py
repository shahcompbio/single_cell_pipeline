'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os

from single_cell.utils import bamutils
from single_cell.utils import helpers


def cell_region_merge_bams(cell_bams, region_bam, region, docker_image):
    cell_bams = cell_bams.values()
    region = '{}:{}-{}'.format(*region.split('-'))

    bamutils.bam_merge(
        cell_bams, region_bam,
        region=region,
        docker_image=docker_image)

    bamutils.bam_index(
        region_bam, region_bam + '.bai',
        docker_image=docker_image)


def merge_bams(bams, outputs, regions, samtools_docker, tempdir, ncores=None):
    merge_tempdir = os.path.join(tempdir, "merge")
    commands = []
    for region in regions:
        output = outputs[region]
        region = '{}:{}-{}'.format(*region.split('-'))
        cmd = list(['samtools', 'merge', '-f', '-R', region])
        cmd.append(output)
        cmd.extend(bams.values())
        commands.append(cmd)
    helpers.run_in_gnu_parallel(commands, merge_tempdir, samtools_docker, ncores=ncores)

    index_tempdir = os.path.join(tempdir, "index")
    commands = []
    for region in regions:
        output = outputs[region]
        commands.append(['samtools', 'index', output, output + ".bai"])

    helpers.run_in_gnu_parallel(commands, index_tempdir, samtools_docker, ncores=ncores)
