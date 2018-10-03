'''
Created on Feb 22, 2018

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
from workflows import merge_bams
from single_cell.utils import helpers

def merge_bams_workflow(workflow, args):


    input_yaml = args["input_yaml"]
    output_template = args["merged_bam_template"]

    info_file = os.path.join(args["out_dir"], "info.yaml")
    config = helpers.load_config(args)
    bam_files, bai_files  = helpers.get_bams(input_yaml)
    cellids = helpers.get_samples(input_yaml)

    wgs_bam_template = output_template
    wgs_bai_template = wgs_bam_template + ".bai"

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cellids,
    )

    workflow.transform(
        name="get_regions",
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 },
        func="single_cell.utils.pysamutils.get_regions_from_reference",
        ret=pypeliner.managed.TempOutputObj('region'),
        args=(
              config["ref_genome"],
              config["split_size"],
              config["chromosomes"],
        )
    )

    workflow.subworkflow(
        name="wgs_merge_workflow",
        func=merge_bams.create_merge_bams_workflow,
        args=(
            mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_files),
            mgd.InputFile('bam_markdups', 'cell_id', fnames=bai_files),
            mgd.OutputFile("merged_bam", "region", axes_origin=[], template=wgs_bam_template),
            mgd.OutputFile("merged_bai", "region", axes_origin=[], template=wgs_bai_template),
            cellids,
            config,
            mgd.TempInputObj("region"),
            mgd.OutputFile(info_file)
        )
    )


    return workflow