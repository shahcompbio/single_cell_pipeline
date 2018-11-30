'''
Created on Jul 11, 2017

@author: dgrewal
'''

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers
import single_cell


def create_merge_bams_workflow(
    input_bams,
    merged_bams,
    cell_ids,
    config,
    regions
):

    baseimage = config['docker']['single_cell_pipeline']

    merged_bams = dict([(region, merged_bams[region])
                         for region in regions])

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cell_ids,
    )

    workflow.setobj(
        obj=mgd.OutputChunks('region'),
        value=regions,
    )

    one_split_job = config["one_split_job"]

    if one_split_job:
        workflow.transform(
            name='merge_bams',
            ctx={'mem': config['memory']['high'], 'ncpus': config['max_cores'], 'docker_image': baseimage},
            func="single_cell.workflows.merge_bams.tasks.merge_bams",
            args=(
                mgd.InputFile('bam', 'cell_id', fnames=input_bams, extensions=['.bai']),
                mgd.OutputFile('merged.bam', "region", fnames=merged_bams, axes_origin=[]),
                regions,
                config['docker']['samtools'],
            ),
            kwargs={"ncores": config["max_cores"]}
        )
    else:
        workflow.transform(
            name='split_merge_tumour',
            func='single_cell.workflows.merge_bams.tasks.cell_region_merge_bams',
            axes=('region',),
            args=(
                mgd.InputFile('tumour_cells.bam', 'cell_id', extensions=['.bai'], fnames=input_bams),
                mgd.OutputFile(
                    'tumour_regions.bam', 'region', axes_origin=[], extensions=['.bai'], fnames=merged_bams),
                mgd.Instance('region'),
                config['docker']['samtools'],
            ),
        )

    return workflow


def create_cell_region_merge_workflow(
    cell_bams,
    region_bams,
    regions,
    docker_image,
):
    region_bams = dict([(region, region_bams[region]) for region in regions])

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cell_bams.keys(),
    )

    workflow.setobj(
        obj=mgd.OutputChunks('region'),
        value=regions,
    )

    workflow.transform(
        name='split_merge_tumour',
        func='single_cell.workflows.merge_bams.tasks.cell_region_merge_bams',
        axes=('region',),
        args=(
            mgd.InputFile('tumour_cells.bam', 'cell_id', extensions=['.bai'], fnames=cell_bams),
            mgd.OutputFile('tumour_regions.bam', 'region', axes_origin=[], extensions=['.bai'], fnames=region_bams),
            mgd.Instance('region'),
            docker_image
        ),
    )

    return workflow

