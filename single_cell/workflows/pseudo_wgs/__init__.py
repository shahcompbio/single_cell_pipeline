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
    sample_ids,
    config,
    regions):
 
 
    merged_bams = dict([(region, merged_bams[region])
                         for region in regions])

    merged_bais = dict([(region, merged_bais[region])
                         for region in regions])


    raise Exception(merged_bais)
 
    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )

    workflow.setobj(
        obj=mgd.OutputChunks('regions'),
        value=regions,
    )


    workflow.transform(
        name='merge_bams',
        ctx={'mem': config["memory"]['high'], 'pool_id': config['pools']['highmem'], 'ncpus':config["max_cores"]},
        func=tasks.merge_bams,
        args=(
            mgd.InputFile('bam', 'sample_id', fnames=input_bams),
            mgd.OutputFile('merged.bam', "regions", fnames=merged_bams, axes_origin=[]),
            mgd.TempOutputFile('merged.bam.bai', "regions", fnames=merged_bais, axes_origin=[]),
            regions
        ),
        kwargs = {"ncores": config["max_cores"]}
    )

    return workflow
