'''
Created on Apr 6, 2018

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
from workflows import split_bams
from single_cell.utils import helpers
import single_cell

def split_bam_workflow(args):

    workflow = pypeliner.workflow.Workflow()
    config = helpers.load_config(args)
    config = config['split_bam']

    baseimage = config['docker']['single_cell_pipeline']

    split_bam_template = args["split_bam_template"]

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
            ctx={'mem': config['memory']['low'], 'ncpus': 1, 'docker_image': baseimage},
            func="single_cell.utils.pysamutils.get_regions_from_reference",
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
            mgd.OutputFile(
                "normal.split.bam", splitkeyword,
                template=split_bam_template, axes_origin=[]
            ),
            pypeliner.managed.TempInputObj(splitkeyword),
            config,
        ),
        kwargs={"by_reads": by_reads}
    )

    return workflow


def split_bam_pipeline(args):

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = split_bam_workflow(args)

    pyp.run(workflow)