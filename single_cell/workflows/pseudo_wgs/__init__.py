'''
Created on Jul 11, 2017

@author: dgrewal
'''


'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd

import tasks


def create_wgs_workflow(
    input_bams,
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
        ctx={'mem': config["memory"]['high'], 'pool_id': config['pools']['highmem'], 'ncpus':config["max_cores"]},
        func=tasks.merge_bams,
        args=(
            mgd.InputFile('bam', 'cell_id', fnames=input_bams),
            mgd.OutputFile('merged.bam', "region", fnames=merged_bams, axes_origin=[]),
            mgd.OutputFile('merged.bam.bai', "region", fnames=merged_bais, axes_origin=[]),
            regions
        ),
        kwargs = {"ncores": config["max_cores"]}
    )

    return workflow
