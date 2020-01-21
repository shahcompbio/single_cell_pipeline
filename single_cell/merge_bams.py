'''
Created on Feb 22, 2018

@author: dgrewal
'''
import os
import sys

import pypeliner
import pypeliner.managed as mgd
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

    bam_files = inpututils.load_merge_cell_bams(args['input_yaml'])

    merge_out_template = os.path.join(args['out_dir'], '{region}.bam')

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=list(bam_files.keys()),
    )

    workflow.transform(
        name="get_regions",
        func="single_cell.utils.pysamutils.get_regions_from_reference",
        ret=pypeliner.managed.OutputChunks('region'),
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
            mgd.InputChunks("region"),
            config,
        )
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='single_cell.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args['out_dir'],
            mgd.Template('bam_filenames', 'region', template=merge_out_template),
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': inpututils.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'template': (mgd.InputChunks('region'), merge_out_template, 'region'),
            'metadata': {
                'type': 'pseudowgs_regionbams',
                'cell_ids': list(bam_files.keys())}

        }
    )

    return workflow


def merge_bams_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = merge_bams_workflow(args)

    pyp.run(workflow)
