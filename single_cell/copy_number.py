"""
Created on Apr 13, 2018

@author: dgrewal
"""

import os
import pypeliner.managed as mgd
from workflows import titan
from single_cell.utils import helpers
from workflows import extract_seqdata
import pypeliner

def copy_number_calling_workflow(args):

    config = helpers.load_config(args)
    config = config['copy_number_calling']

    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1,
           'docker_image': config['docker']['single_cell_pipeline']
    }
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    data = helpers.load_pseudowgs_input(args['input_yaml'])
    normal_wgs = data['normal_wgs']
    tumour_cells = data['tumour_cells']
    assert '{region}' in normal_wgs


    copynumber_dir = os.path.join(args["out_dir"], "copynumber")

    out_file = os.path.join(copynumber_dir, "results", "results.h5")

    cloneid = args["clone_id"]

    remixt_config = config.get('extract_seqdata', {})

    workflow.setobj(
        obj=mgd.OutputChunks('tumour_cell_id'),
        value=tumour_cells.keys(),
    )

    workflow.transform(
        name="get_regions",
        ctx=dict(mem=config['memory']['low']),
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
        func="remixt.config.get_param",
        ret=mgd.TempOutputObj('bam_max_fragment_length'),
        args=(
              remixt_config,
              'bam_max_fragment_length'
        )
    )

    workflow.transform(
        name="get_bam_max_soft_clipped",
        func="remixt.config.get_param",
        ret=mgd.TempOutputObj('bam_max_soft_clipped'),
        args=(
              remixt_config,
              'bam_max_soft_clipped'
        )
    )

    workflow.transform(
        name="get_bam_check_proper_pair",
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
                fnames=tumour_cells,
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
        ctx={'disk': 200},
        func=extract_seqdata.create_extract_seqdata_workflow,
        args=(
            mgd.InputFile(
                'bam_markdups',
                'region',
                template=normal_wgs,
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
            tumour_cells.keys(),
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
