'''
Created on Apr 6, 2018

@author: dgrewal
'''
import os

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import inpututils
from single_cell.utils import helpers
from single_cell.workflows import split_bams


def split_bam_workflow(args):
    workflow = pypeliner.workflow.Workflow()
    config = inpututils.load_config(args)
    config = config['split_bam']

    _, _, bam_file = inpututils.load_split_wgs_input(args['input_yaml'])

    baseimage = config['docker']['single_cell_pipeline']

    split_bam_template = os.path.join(args['out_dir'], '{region}_split_wgs.bam')

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
        ctx={'mem': config['memory']['low'], 'ncpus': 1, 'docker_image': baseimage},
        args=(
            mgd.InputFile(bam_file),
            mgd.OutputFile(
                "normal.split.bam", 'region',
                template=split_bam_template, axes_origin=[]
            ),
            pypeliner.managed.TempInputObj('region'),
            config,
        ),
    )

    return workflow


def generate_meta_files(args):
    out_dir = args["out_dir"]

    sample_id, library_id, bam_files = inpututils.load_split_wgs_input(args['input_yaml'])

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



def split_bam_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = split_bam_workflow(args)

    pyp.run(workflow)
