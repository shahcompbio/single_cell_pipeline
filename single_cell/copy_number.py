'''
Created on Apr 13, 2018

@author: dgrewal
'''

import os
import pypeliner.managed as mgd
from workflows import titan
from single_cell.utils import helpers
from workflows import extract_seqdata


def copy_number_calling_workflow(workflow, args):

    config = helpers.load_config(args)

    ctx = {'mem_retry_increment': 2, 'ncpus': 1,
           'mem': config["memory"]['low'],
           'pool_id': config['pools']['standard']}
    docker_ctx = helpers.get_container_ctx(config['containers'], 'single_cell_pipeline')
    ctx.update(docker_ctx)

    tumour_bam_files, tumour_bai_files = helpers.get_bams(args['tumour_yaml'])

    normal_bam_files, normal_bai_files = helpers.get_bams(args['normal_yaml'])

    tumour_cellids = helpers.get_samples(args['tumour_yaml'])

    normal_cellids = helpers.get_samples(args['normal_yaml'])

    if set(tumour_bam_files.keys()) != set(tumour_cellids):
        raise ValueError()

    if set(normal_bam_files.keys()) != set(normal_cellids):
        raise ValueError()

    copynumber_dir = os.path.join(args["out_dir"], "copynumber")

    out_file = os.path.join(copynumber_dir, "results", "results.h5")

    cloneid = args["clone_id"]

    remixt_config = config['titan_params'].get('extract_seqdata', {})

    workflow.setobj(
        obj=mgd.OutputChunks('tumour_cell_id'),
        value=tumour_cellids,
    )

    workflow.setobj(
        obj=mgd.OutputChunks('normal_cell_id'),
        value=normal_cellids,
    )

    workflow.transform(
        name="get_snp_positions_filename",
        ctx=ctx,
        func="remixt.config.get_filename",
        ret=mgd.TempOutputObj('snp_positions_filename'),
        args=(
              remixt_config,
              config['titan_params']['ref_data_dir'],
              'snp_positions'
        )
    )

    workflow.transform(
        name="get_bam_max_fragment_length",
        ctx=ctx,
        func="remixt.config.get_param",
        ret=mgd.TempOutputObj('bam_max_fragment_length'),
        args=(
              remixt_config,
              'bam_max_fragment_length'
        )
    )

    workflow.transform(
        name="get_bam_max_soft_clipped",
        ctx=ctx,
        func="remixt.config.get_param",
        ret=mgd.TempOutputObj('bam_max_soft_clipped'),
        args=(
              remixt_config,
              'bam_max_soft_clipped'
        )
    )

    workflow.transform(
        name="get_bam_check_proper_pair",
        ctx=ctx,
        func="remixt.config.get_param",
        ret=mgd.TempOutputObj('bam_check_proper_pair'),
        args=(
              remixt_config,
              'bam_check_proper_pair'
        )
    )


    workflow.subworkflow(
        name="extract_seqdata_tumour",
        axes=('tumour_cell_id',),
        func=extract_seqdata.create_extract_seqdata_workflow,
        args=(
            mgd.InputFile(
                'bam_markdups',
                'tumour_cell_id',
                fnames=tumour_bam_files),
            mgd.InputFile(
                'bam_markdups_index',
                'tumour_cell_id',
                fnames=tumour_bai_files),
            mgd.TempOutputFile("tumour.h5", "tumour_cell_id"),
            config,
            config['titan_params'].get('extract_seqdata', {}),
            config['titan_params']['ref_data_dir'],
            mgd.TempInputObj('snp_positions_filename'),
            mgd.TempInputObj('bam_max_fragment_length'),
            mgd.TempInputObj('bam_max_soft_clipped'),
            mgd.TempInputObj('bam_check_proper_pair'),
        )
    )

    workflow.subworkflow(
        name="extract_seqdata_normal",
        axes=('normal_cell_id',),
        func=extract_seqdata.create_extract_seqdata_workflow,
        args=(
            mgd.InputFile(
                'bam_markdups',
                'normal_cell_id',
                fnames=normal_bam_files),
            mgd.InputFile(
                'bam_markdups_index',
                'normal_cell_id',
                fnames=normal_bai_files),
            mgd.TempOutputFile("normal.h5", "normal_cell_id"),
            config,
            config['titan_params'].get('extract_seqdata', {}),
            config['titan_params']['ref_data_dir'],
            mgd.TempInputObj('snp_positions_filename'),
            mgd.TempInputObj('bam_max_fragment_length'),
            mgd.TempInputObj('bam_max_soft_clipped'),
            mgd.TempInputObj('bam_check_proper_pair'),
        )
    )

    workflow.subworkflow(
        name='titan_workflow',
        func=titan.create_titan_workflow,
        args=(
            mgd.TempInputFile("normal.h5", "normal_cell_id"),
            mgd.TempInputFile("tumour.h5", "tumour_cell_id"),
            config['ref_genome'],
            copynumber_dir,
            out_file,
            config,
            args,
            tumour_cellids,
            normal_cellids,
            cloneid
        ),
    )

    return workflow
