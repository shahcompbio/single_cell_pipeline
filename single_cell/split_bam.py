'''
Created on Apr 6, 2018

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
from workflows import split_bams
from single_cell.utils import helpers

def split_bam_workflow(workflow, args):

    args = helpers.generate_configs_in_temp(args)

    config = helpers.load_config(args)

    normal_bam_template = args["split_bam_template"]
    normal_bai_template = args["split_bam_template"] + ".bai"


    workflow.transform(
        name="get_regions",
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 },
        func=helpers.get_bam_regions,
        ret=pypeliner.managed.TempOutputObj('region'),
        args=(
              config["ref_genome"],
              config["split_size"],
              config["chromosomes"],
        )
    )

    workflow.subworkflow(
        name="split_normal",
        func=split_bams.create_split_workflow,
        args = (
            mgd.InputFile(args['wgs_bam']),
            mgd.InputFile(args['wgs_bam'] + ".bai"),
            mgd.OutputFile("normal.split.bam", "region", template=normal_bam_template, axes_origin=[]),
            mgd.OutputFile("normal.split.bam.bai", "region", template=normal_bai_template, axes_origin=[]),
            pypeliner.managed.TempInputObj('region'),
            config,
        )
    )

    return workflow