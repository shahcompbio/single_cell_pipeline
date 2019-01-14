"""
Created on Apr 13, 2018

@author: dgrewal
"""

import os
import pypeliner.managed as mgd
from workflows import titan
from single_cell.utils import helpers
from workflows import extract_seqdata

def copy_number_calling_workflow(workflow, args):

    config = helpers.load_config(args)
    config = config['copy_number_calling']

    ctx = {'mem_retry_increment': 2, 'ncpus': 1,
           'docker_image': config['docker']['single_cell_pipeline']
    }

    tumour_bam_files, tumour_bai_files = helpers.get_bams(args['tumour_yaml'])

    tumour_cellids = helpers.get_samples(args['tumour_yaml'])

    if set(tumour_bam_files.keys()) != set(tumour_cellids):
        raise ValueError()

    copynumber_dir = os.path.join(args["out_dir"], "copynumber")

    out_file = os.path.join(copynumber_dir, "results", "results.h5")

    cloneid = args["clone_id"]

    remixt_config = config.get('extract_seqdata', {})

    workflow.setobj(
        obj=mgd.OutputChunks('tumour_cell_id'),
        value=tumour_cellids,
    )

    workflow.transform(
        name="get_regions",
        ctx=dict(mem=config['memory']['low'], **ctx),
        func="single_cell.utils.pysamutils.get_regions_from_reference",
        ret=mgd.OutputChunks('region'),
        args=(
            config["ref_genome"],
            config["split_size"],
            config["chromosomes"],
        )
    )

    workflow.transform(
        name="get_snp_positions_filename",
        ctx=ctx,
        func="remixt.config.get_filename",
        ret=mgd.TempOutputObj('snp_positions_filename'),
        args=(
              remixt_config,
              config['ref_data_dir'],
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
                fnames=tumour_bam_files,
                extensions=['.bai']
            ),
            mgd.TempOutputFile("tumour.h5", "tumour_cell_id"),
            config.get('extract_seqdata', {}),
            config['ref_data_dir'],
            config
        )
    )

    workflow.subworkflow(
        name="extract_seqdata_normal",
        axes=('region',),
        func=extract_seqdata.create_extract_seqdata_workflow,
        args=(
            mgd.InputFile(
                'bam_markdups',
                'region',
                template=args['normal_template'],
                extensions=['.bai']
            ),
            mgd.TempOutputFile("normal.h5", "region"),
            config.get('extract_seqdata', {}),
            config['ref_data_dir'],
            config,
        )
    )

    workflow.subworkflow(
        name='titan_workflow',
        func=titan.create_titan_workflow,
        args=(
            mgd.TempInputFile("normal.h5", "region"),
            mgd.TempInputFile("tumour.h5", "tumour_cell_id"),
            config['ref_genome'],
            copynumber_dir,
            mgd.OutputFile(out_file),
            config,
            args,
            tumour_cellids,
            mgd.InputChunks('region'),
            cloneid
        ),
    )
    #
    # info_file = os.path.join(args["out_dir"],'results','copynumber_calling', "info.yaml")
    #
    # results = {
    #     'copynumber_data': helpers.format_file_yaml(out_file),
    # }
    #
    # tumours = {k: helpers.format_file_yaml(v) for k,v in tumour_bam_files.iteritems()}
    # normals = {k: helpers.format_file_yaml(v) for k,v in normal_bam_files.iteritems()}
    # input_datasets = {'tumour': tumours, 'normal': normals}
    #
    # metadata = {
    #     'copynumber_calling': {
    #         'chromosomes': config['chromosomes'],
    #         'ref_genome': config['ref_genome'],
    #         'version': single_cell.__version__,
    #         'results': results,
    #         'containers': config['containers'],
    #         'input_datasets': input_datasets,
    #         'output_datasets': None
    #     }
    # }
    #
    # workflow.transform(
    #     name='generate_meta_yaml',
    #     ctx=dict(mem=config['memory']['med'],
    #              pool_id=config['pools']['standard'],
    #              mem_retry_increment=2, ncpus=1),
    #     func="single_cell.utils.helpers.write_to_yaml",
    #     args=(
    #         mgd.OutputFile(info_file),
    #         metadata
    #     )
    # )

    return workflow
