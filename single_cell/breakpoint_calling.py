import os

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import inpututils


def breakpoint_calling_workflow(args):
    run_destruct = args['destruct']
    run_lumpy = args['lumpy']
    if not any((run_destruct, run_lumpy)):
        run_destruct = True
        run_lumpy = True

    config = inpututils.load_config(args)
    config = config['breakpoint_calling']

    normal_data, tumour_cells = inpututils.load_breakpoint_calling_input(args['input_yaml'])

    bkp_dir = os.path.join(args['out_dir'], 'breakpoint_calling')

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
        breakpoints_filename = os.path.join(bkp_dir, 'breakpoints.csv.gz')
        breakpoints_lib_filename = os.path.join(bkp_dir, 'breakpoints_library.csv.gz')
        cell_counts_filename = os.path.join(bkp_dir, 'cell_counts.csv.gz')

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
                mgd.OutputFile(breakpoints_filename, extensions=['.yaml']),
                mgd.OutputFile(breakpoints_lib_filename, extensions=['.yaml']),
                mgd.OutputFile(cell_counts_filename, extensions=['.yaml']),
            ),
        )

    if run_lumpy:
        breakpoints_bed = os.path.join(bkp_dir, 'lumpy_breakpoints.bed')
        breakpoints_csv = os.path.join(bkp_dir, 'lumpy_breakpoints.csv.gz')
        breakpoints_evidence_csv = os.path.join(bkp_dir, 'lumpy_breakpoints_evidence.csv.gz')

        workflow.subworkflow(
            name='lumpy',
            func="single_cell.workflows.lumpy.create_lumpy_workflow",
            args=(
                config,
                normal_bam,
                mgd.InputFile('tumour.bam', 'tumour_cell_id', fnames=tumour_cells, extensions=['.bai']),
                mgd.OutputFile(breakpoints_csv, extensions=['.yaml']),
                mgd.OutputFile(breakpoints_evidence_csv, extensions=['.yaml']),
                mgd.OutputFile(breakpoints_bed),
            ),
        )

    return workflow


def breakpoint_calling_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = breakpoint_calling_workflow(args)

    pyp.run(workflow)
