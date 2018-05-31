'''
Created on Jul 11, 2017

@author: dgrewal
'''

import pypeliner
import pypeliner.managed as mgd


def create_wgs_workflow(
    input_bams,
    input_bais,
    merged_bams,
    merged_bais,
    cell_ids,
    config,
    regions):
 
 
    merged_bams = dict([(region, merged_bams[region])
                         for region in regions])

    merged_bais = dict([(region, merged_bais[region])
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


    workflow.transform(
        name='merge_bams',
        ctx={'mem': config["memory"]['high'], 'pool_id': config['pools']['multicore'], 'ncpus':config["max_cores"]},
        func="single_cell.workflows.pseudo_wgs.tasks.merge_bams",
        args=(
            mgd.InputFile('bam', 'cell_id', fnames=input_bams),
            mgd.InputFile('bai', 'cell_id', fnames=input_bais),
            mgd.OutputFile('merged.bam', "region", fnames=merged_bams, axes_origin=[]),
            mgd.OutputFile('merged.bam.bai', "region", fnames=merged_bais, axes_origin=[]),
            regions
        ),
        kwargs = {"ncores": config["max_cores"]}
    )

    return workflow
