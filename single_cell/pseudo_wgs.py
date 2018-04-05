'''
Created on Feb 22, 2018

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
from workflows import pseudo_wgs 
from single_cell.utils import helpers

def pseudo_wgs_workflow(workflow, args):

    config = helpers.load_config(args)
    bam_files, _  = helpers.get_bams(args['input_yaml'])
    cellids = helpers.get_samples(args['input_yaml'])

    wgs_bam_template = args["merged_wgs_template"]
    wgs_bai_template = args["merged_wgs_template"] + ".bai"



    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cellids,
    )

    workflow.transform(
        name="get_regions",
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 },
        func=helpers.get_bam_regions,
        ret=pypeliner.managed.TempOutputObj('regions'),
        args=(
              config["ref_genome"],
              config["split_size"],
              config["chromosomes"],
        )
    )

    workflow.subworkflow(
        name='wgs_workflow',
        func=pseudo_wgs.create_wgs_workflow,
        args=(
            mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_files),
            mgd.OutputFile("merged_bam", "regions", axes_origin=[], template=wgs_bam_template),
            mgd.OutputFile("merged_bai", "regions", axes_origin=[], template=wgs_bai_template),
            cellids,
            config,
            mgd.TempInputObj("regions"),
        )
    )


    return workflow