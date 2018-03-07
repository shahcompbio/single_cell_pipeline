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
    sampleids = helpers.get_samples(args['input_yaml'])

    pseudo_wgs_bam = args["merged_wgs"]
    pseudo_wgs_bai = args["merged_wgs"]+".bai"


    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sampleids,
    )
    
    workflow.subworkflow(
        name='wgs_workflow',
        func=pseudo_wgs.create_wgs_workflow,
        args=(
            mgd.InputFile('bam_markdups', 'sample_id', fnames=bam_files),
            mgd.OutputFile(pseudo_wgs_bam),
            mgd.OutputFile(pseudo_wgs_bai),
            config['ref_genome'],
            sampleids,
            config,
            args['out_dir'],
        )
    )


    return workflow