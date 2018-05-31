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

    config = helpers.load_config(args)

    split_bam_template = args["split_bam_template"]
    split_bai_template = args["split_bam_template"] + ".bai"

    by_reads = False if "{region}" in split_bam_template else True
    splitkeyword = "region" if "{region}" in split_bam_template else "reads"

    if by_reads:
        splitnames = [str(i) for i in range(config["num_splits_byreads"])]

        workflow.setobj(
            obj=mgd.OutputChunks('reads'),
            value=splitnames,
        )

    else:
        workflow.transform(
            name="get_regions",
            ctx={
                'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2,
                'pool_id': config['pools']['standard'], 'ncpus': 1
            },
            func=helpers.get_regions_from_reference,
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
        args=(
            mgd.InputFile(args['wgs_bam']),
            mgd.InputFile(args['wgs_bam'] + ".bai"),
            mgd.OutputFile(
                "normal.split.bam", splitkeyword,
                template=split_bam_template, axes_origin=[]
            ),
            mgd.OutputFile(
                "normal.split.bam.bai", splitkeyword,
                template=split_bai_template, axes_origin=[]
            ),
            pypeliner.managed.TempInputObj(splitkeyword),
            config,
        ),
        kwargs={"by_reads": by_reads}
    )

    return workflow
