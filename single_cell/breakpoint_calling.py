import os
import sys

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import inpututils


def get_output_files(bkp_dir, run_destruct, run_lumpy):
    data = {}

    if run_destruct:
        data['destruct_breakpoints_filename'] = os.path.join(bkp_dir, 'destruct_breakpoints.csv.gz')
        data['destruct_breakpoints_lib_filename'] = os.path.join(bkp_dir, 'destruct_breakpoints_library.csv.gz')
        data['destruct_cell_counts_filename'] = os.path.join(bkp_dir, 'destruct_cell_counts.csv.gz')

    if run_lumpy:
        data['lumpy_breakpoints_bed'] = os.path.join(bkp_dir, 'lumpy_breakpoints.bed')
        data['lumpy_breakpoints_csv'] = os.path.join(bkp_dir, 'lumpy_breakpoints.csv.gz')
        data['lumpy_breakpoints_evidence_csv'] = os.path.join(bkp_dir, 'lumpy_breakpoints_evidence.csv.gz')

    return data


def breakpoint_calling_workflow(args):
    config = inpututils.load_config(args)
    config = config['breakpoint_calling']

    run_destruct = True if args['destruct'] else False
    run_lumpy = True if args['lumpy'] else False

    if not run_destruct and not run_lumpy:
        run_destruct = True
        run_lumpy = True

    normal_data, tumour_cells = inpututils.load_breakpoint_calling_input(args['input_yaml'])

    bkp_dir = os.path.join(args['out_dir'])
    bkp_meta = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    out_files = get_output_files(bkp_dir, run_destruct, run_lumpy)

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

    if run_destruct:
        workflow.subworkflow(
            name='destruct',
            ctx={'docker_image': config['docker']['single_cell_pipeline']},
            func="single_cell.workflows.destruct_singlecell.create_destruct_workflow",
            args=(
                normal_bam,
                mgd.InputFile('tumour.bam', 'tumour_cell_id', fnames=tumour_cells),
                config.get('destruct_config', {}),
                config,
                ref_data_directory,
                mgd.OutputFile(out_files['destruct_breakpoints_filename'], extensions=['.yaml']),
                mgd.OutputFile(out_files['destruct_breakpoints_lib_filename'], extensions=['.yaml']),
                mgd.OutputFile(out_files['destruct_cell_counts_filename'], extensions=['.yaml']),
            ),
        )

    if run_lumpy:
        workflow.subworkflow(
            name='lumpy',
            func="single_cell.workflows.lumpy.create_lumpy_workflow",
            args=(
                config,
                normal_bam,
                mgd.InputFile('tumour.bam', 'tumour_cell_id', fnames=tumour_cells, extensions=['.bai']),
                mgd.OutputFile(out_files['lumpy_breakpoints_csv'], extensions=['.yaml']),
                mgd.OutputFile(out_files['lumpy_breakpoints_evidence_csv'], extensions=['.yaml']),
                mgd.OutputFile(out_files['lumpy_breakpoints_bed']),
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
