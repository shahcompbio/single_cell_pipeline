import os
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers
import single_cell

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


    info_file = os.path.join(args["out_dir"],'results','breakpoint_calling', "info.yaml")

    results = {
        'destruct_data': helpers.format_file_yaml(breakpoints_filename),
    }

    input_datasets = {k: helpers.format_file_yaml(v) for k,v in bam_files.iteritems()}
    input_datasets = {'normal': normal_bam_file,
                      'tumour': input_datasets}

    metadata = {
        'breakpoint_calling': {
            'ref_data': ref_data_directory,
            'version': single_cell.__version__,
            'results': results,
            'containers': config['containers'],
            'input_datasets': input_datasets,
            'output_datasets': None
        }
    }

    workflow.transform(
        name='generate_meta_yaml',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 mem_retry_increment=2, ncpus=1),
        func="single_cell.utils.helpers.write_to_yaml",
        args=(
            mgd.OutputFile(info_file),
            metadata
        )
    )

    return workflow

