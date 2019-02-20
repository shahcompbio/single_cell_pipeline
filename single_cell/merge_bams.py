'''
Created on Feb 22, 2018

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
from workflows import merge_bams
from single_cell.utils import helpers
import single_cell

def merge_bams_workflow(workflow, args):

    config = helpers.load_config(args)
    config = config['merge_bams']

    baseimage = config['docker']['single_cell_pipeline']

    data = helpers.load_pseudowgs_input(args['input_yaml'])
    tumour_wgs = data['tumour_wgs']
    normal_wgs = data['normal_wgs']
    tumour_cells = data['tumour_cells']
    normal_cells = data['normal_cells']

    bam_files = tumour_cells if tumour_cells else normal_cells
    wgs_bams = tumour_wgs if tumour_cells else normal_wgs


    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=bam_files.keys(),
    )

    if isinstance(wgs_bams, dict):
        workflow.setobj(
            obj=mgd.OutputChunks('regions'),
            value=wgs_bams.keys(),
        )
        workflow.set_filenames("merged.bam", "region", fnames=wgs_bams)
    else:
        workflow.transform(
            name="get_regions",
            ctx={'mem_retry_increment': 2, 'ncpus': 1, 'mem': config["memory"]['low'], 'docker_image': baseimage},
            func="single_cell.utils.pysamutils.get_regions_from_reference",
            ret=pypeliner.managed.OutputChunks('region'),
            args=(
                config["ref_genome"],
                config["split_size"],
                config["chromosomes"],
            )
        )
        workflow.set_filenames('merged.bam', 'region', template=wgs_bams)

    workflow.subworkflow(
        name="wgs_merge_workflow",
        func=merge_bams.create_merge_bams_workflow,
        args=(
            mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_files, extensions=['.bai']),
            mgd.OutputFile("merged.bam", "region", axes_origin=[], extensions=['.bai']),
            mgd.TempInputObj("region"),
            config,
        )
    )

    workflow.transform(
        name="get_files",
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func='single_cell.utils.helpers.resolve_template',
        ret=pypeliner.managed.TempOutputObj('outputs'),
        args=(
            pypeliner.managed.TempInputObj('region'),
            wgs_bams,
            'region'
        )

    )

    info_file = os.path.join(args["out_dir"], 'results','merge_bams', "info.yaml")

    inputs = {k: helpers.format_file_yaml(v) for k,v in bam_files.iteritems()}

    metadata = {
        'merge_bams': {
            'name': 'merge_bams',
            'ref_genome': config["ref_genome"],
            'version': single_cell.__version__,
            'containers': config['docker'],
            'output_datasets': pypeliner.managed.TempInputObj('outputs'),
            'input_datasets': inputs,
            'results': None
        }
    }

    workflow.transform(
        name='generate_meta_yaml',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.utils.helpers.write_to_yaml",
        args=(
            mgd.OutputFile(info_file),
            metadata
        )
    )

    return workflow


