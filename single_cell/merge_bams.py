'''
Created on Feb 22, 2018

@author: dgrewal
'''
import os

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers
from single_cell.utils import inpututils
from single_cell.workflows import merge_bams


def merge_bams_workflow(args):
    config = inpututils.load_config(args)
    config = config['merge_bams']

    baseimage = config['docker']['single_cell_pipeline']

    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50,
           'ncpus': 1, 'mem': config["memory"]['low'],
           'docker_image': baseimage}
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    _, _, bam_files = inpututils.load_merge_cell_bams(args['input_yaml'])

    merge_out_template = os.path.join(args['out_dir'], '{region}.bam')

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=list(bam_files.keys()),
    )

    workflow.transform(
        name="get_regions",
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
            mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_files, extensions=['.bai']),
            mgd.OutputFile("merged.bam", "region", axes_origin=[], extensions=['.bai'], template=merge_out_template),
            mgd.TempInputObj("region"),
            config,
        )
    )

    return workflow


def generate_meta_files(args):
    out_dir = args["out_dir"]

    sample_id, library_id, bam_files = inpututils.load_merge_cell_bams(args['input_yaml'])

    metadata = dict(
        library_id=library_id, sample_ids=sample_id,
    )

    # filepaths = helpers.resolve_template(os.path.join(args['out_dir'], '{region}.bam'), regions)
    filepaths = []

    metadata['type'] = 'merge_cell_bams'
    helpers.generate_and_upload_metadata(
        args,
        out_dir,
        [],
        metadata,
        input_yaml=args['input_yaml']
    )


def merge_bams_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = merge_bams_workflow(args)

    pyp.run(workflow)

    generate_meta_files(args)
