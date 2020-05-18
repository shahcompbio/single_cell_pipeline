'''
Created on Aug 29, 2018

@author: dgrewal
'''

import os
import sys

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import inpututils


def infer_haps_workflow(args):
    config = inpututils.load_config(args)
    config = config['infer_haps']
    baseimage = config['docker']['single_cell_pipeline']

    ctx = dict(mem_retry_increment=2, disk_retry_increment=50, ncpus=1, docker_image=baseimage)
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    haplotypes_filename = os.path.join(args["out_dir"], "haplotypes.csv.gz")

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    normal_data = inpututils.load_infer_haps_input(args['input_yaml'])

    if isinstance(normal_data, dict):
        workflow.setobj(
            obj=mgd.OutputChunks('normal_cell_id'),
            value=list(normal_data.keys()),
        )
        bam_file = mgd.InputFile('normal.bam', 'normal_cell_id', fnames=normal_data, extensions=['.bai'])
    else:
        bam_file = mgd.InputFile(normal_data, extensions=['.bai'])

    workflow.subworkflow(
        name='infer_haps',
        func='single_cell.workflows.infer_haps.infer_haps',
        args=(
            bam_file,
            mgd.OutputFile(haplotypes_filename, extensions=['.yaml']),
            config,
        ),
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='single_cell.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args['out_dir'],
            [haplotypes_filename],
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': inpututils.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'infer_haps'}
        }
    )

    return workflow


def count_haps_workflow(args):
    config = inpututils.load_config(args)
    config = config['count_haps']
    baseimage = config['docker']['single_cell_pipeline']

    ctx = dict(mem_retry_increment=2, disk_retry_increment=50, ncpus=1, docker_image=baseimage)
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    allele_counts_filename = os.path.join(args["out_dir"], "allele_counts.csv.gz")

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    haplotypes_filename, tumour_cells = inpututils.load_count_haps_input(args['input_yaml'])

    workflow.setobj(
        obj=mgd.OutputChunks('tumour_cell_id'),
        value=list(tumour_cells.keys()),
    )

    workflow.subworkflow(
        name='extract_allele_readcounts',
        func='single_cell.workflows.extract_allele_readcounts.extract_allele_readcounts',
        args=(
            mgd.InputFile(haplotypes_filename, extensions=['.yaml']),
            mgd.InputFile('tumour_cells.bam', 'tumour_cell_id', extensions=['.bai'],
                          axes_origin=[], fnames=tumour_cells),
            mgd.OutputFile(allele_counts_filename),
            config,
        ),
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='single_cell.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args['out_dir'],
            [allele_counts_filename],
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': inpututils.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'count_haps'}
        }
    )

    return workflow


def infer_haps_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = infer_haps_workflow(args)

    pyp.run(workflow)


def count_haps_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = count_haps_workflow(args)

    pyp.run(workflow)
