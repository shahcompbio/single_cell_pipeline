import os
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers


def breakpoint_calling_workflow(workflow, args):

    config = helpers.load_config(args)

    normal_bam_file = args['matched_normal']
    bam_files, bai_files  = helpers.get_bams(args['input_yaml'])

    varcalls_dir = os.path.join(
        args['out_dir'], 'results', 'breakpoint_calling')
    raw_data_directory = os.path.join(varcalls_dir, 'raw')
    breakpoints_filename = os.path.join(varcalls_dir, 'breakpoints.h5')
    ref_data_directory = '/refdata'

    pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=bam_files.keys(),
    )

    workflow.subworkflow(
        name='destruct',
        func="biowrappers.components.breakpoint_calling.destruct.destruct_pipeline",
        args=(
            mgd.InputFile(normal_bam_file),
            mgd.InputFile('tumour.bam', 'cell_id', fnames=bam_files),
            config.get('destruct', {}),
            ref_data_directory,
            mgd.OutputFile(breakpoints_filename),
            raw_data_directory,
        ),
    )

    return workflow

