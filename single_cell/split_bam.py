'''
Created on Apr 6, 2018

@author: dgrewal
'''
import os
import sys

import pypeliner.managed as mgd
from single_cell.utils import inpututils
from single_cell.workflows import split_bams

import pypeliner


def split_bam_workflow(args):
    config = inpututils.load_config(args)
    config = config['split_bam']

    bam_file = inpututils.load_split_wgs_input(args['input_yaml'])

    baseimage = config['docker']['single_cell_pipeline']

    split_bam_template = os.path.join(args['out_dir'], '{region}.bam')

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    workflow = pypeliner.workflow.Workflow(ctx={'docker_image': baseimage})

    workflow.transform(
        name="get_regions",
        ctx={'mem': config['memory']['low'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.utils.pysamutils.get_regions_from_reference",
        ret=pypeliner.managed.OutputChunks('region'),
        args=(
            config["ref_genome"],
            config["split_size"],
            config["chromosomes"],
        )
    )

    workflow.subworkflow(
        name="split_normal",
        func=split_bams.create_split_workflow,
        ctx={'mem': config['memory']['low'], 'ncpus': 1},
        args=(
            mgd.InputFile(bam_file),
            mgd.OutputFile(
                "normal.split.bam", 'region',
                template=split_bam_template, axes_origin=[]
            ),
            pypeliner.managed.InputChunks('region'),
            config,
        ),
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='single_cell.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args['out_dir'],
            mgd.Template('bam_filenames', 'region', template=split_bam_template),
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': inpututils.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'wgs_regionbams'},
            'template': (mgd.InputChunks('region'), split_bam_template, 'region'),
        }
    )

    return workflow


def split_bam_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = split_bam_workflow(args)

    pyp.run(workflow)
