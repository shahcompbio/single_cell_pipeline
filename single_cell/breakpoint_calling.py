import os
import sys

import pypeliner.managed as mgd
from single_cell.utils import inpututils

import pypeliner


def get_output_files(bkp_dir):
    data = {
        'breakpoints_filename': os.path.join(bkp_dir, 'breakpoints.csv.gz'),
        'breakpoints_lib_filename': os.path.join(bkp_dir, 'breakpoints_library.csv.gz'),
        'cell_counts_filename': os.path.join(bkp_dir, 'cell_counts.csv.gz'),
        'breakpoints_bed': os.path.join(bkp_dir, 'lumpy_breakpoints.bed'),
        'breakpoints_csv': os.path.join(bkp_dir, 'lumpy_breakpoints.csv.gz'),
        'breakpoints_evidence_csv': os.path.join(bkp_dir, 'lumpy_breakpoints_evidence.csv.gz'),

    }

    return data


def breakpoint_calling_workflow(args):
    config = inpututils.load_config(args)
    config = config['breakpoint_calling']

    normal_data, tumour_cells = inpututils.load_breakpoint_calling_input(args['input_yaml'])

    bkp_dir = os.path.join(args['out_dir'])
    bkp_meta = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    out_files = get_output_files(bkp_dir)

    ref_data_directory = config['ref_data_directory']

    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config['docker']['single_cell_pipeline']}
    )

    workflow.setobj(
        obj=mgd.OutputChunks('tumour_cell_id'),
        value=list(tumour_cells.keys()),
    )

    if isinstance(normal_data, dict):
        workflow.setobj(
            obj=mgd.OutputChunks('normal_cell_id'),
            value=list(normal_data.keys()),
        )
        normal_bam = mgd.InputFile(
            'normal_cells.bam', 'normal_cell_id',
            extensions=['.bai'], fnames=normal_data
        )
    else:
        normal_bam = mgd.InputFile(normal_data, extensions=['.bai'])

    workflow.subworkflow(
        name='destruct',
        ctx={'docker_image': config['docker']['destruct']},
        func="single_cell.workflows.destruct_singlecell.create_destruct_workflow",
        args=(
            normal_bam,
            mgd.InputFile('tumour.bam', 'tumour_cell_id', fnames=tumour_cells),
            config.get('destruct_config', {}),
            config,
            ref_data_directory,
            mgd.OutputFile(out_files['breakpoints_filename'], extensions=['.yaml']),
            mgd.OutputFile(out_files['breakpoints_lib_filename'], extensions=['.yaml']),
            mgd.OutputFile(out_files['cell_counts_filename'], extensions=['.yaml']),
        ),
    )

    workflow.subworkflow(
        name='lumpy',
        func="single_cell.workflows.lumpy.create_lumpy_workflow",
        args=(
            config,
            normal_bam,
            mgd.InputFile('tumour.bam', 'tumour_cell_id', fnames=tumour_cells, extensions=['.bai']),
            mgd.OutputFile(out_files['breakpoints_csv'], extensions=['.yaml']),
            mgd.OutputFile(out_files['breakpoints_evidence_csv'], extensions=['.yaml']),
            mgd.OutputFile(out_files['breakpoints_bed']),
        ),
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='single_cell.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            bkp_dir,
            list(out_files.values()),
            mgd.OutputFile(bkp_meta)
        ),
        kwargs={
            'input_yaml_data': inpututils.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'breakpoint_calling'}
        }
    )

    return workflow


def breakpoint_calling_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = breakpoint_calling_workflow(args)

    pyp.run(workflow)
