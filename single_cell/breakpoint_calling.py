import os
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers
import destruct.workflow


default_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']


def breakpoint_calling_workflow(workflow, args):

    config = helpers.load_config(args)

    bam_files, bai_files  = helpers.get_bams(args['input_yaml'])

    varcalls_dir = os.path.join(
        args['out_dir'], 'results', 'breakpoint_calling')
    raw_data_directory = os.path.join(varcalls_dir, 'raw')
    breakpoints_filename = os.path.join(varcalls_dir, 'breakpoints.tsv')
    breakpoint_library_filename = os.path.join(varcalls_dir, 'breakpoint_library.tsv')
    breakpoint_read_filename = os.path.join(varcalls_dir, 'breakpoint_read.tsv')
    ref_data_directory = '/refdata'

    pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=bam_files.keys(),
    )

    workflow.subworkflow(
        name='destruct',
        func=destruct.workflow.create_destruct_workflow,
        args=(
            mgd.InputFile('tumour.bam', 'cell_id', fnames=bam_files),
            mgd.OutputFile(breakpoints_filename),
            mgd.OutputFile(breakpoint_library_filename),
            mgd.OutputFile(breakpoint_read_filename),
            config.get('destruct', {}),
            ref_data_directory,
            raw_data_directory,
        ),
    )

    return workflow

