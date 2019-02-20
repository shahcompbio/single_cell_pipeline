import os
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers
import single_cell


def breakpoint_calling_workflow(workflow, args):

    config = helpers.load_config(args)
    config = config['breakpoint_calling']

    baseimage = config['docker']['single_cell_pipeline']

    tumour_id, tumour_bams, normal_id, normal_bams = helpers.load_pseudowgs_input(args['input_yaml'])

    varcalls_dir = os.path.join(
        args['out_dir'], 'results', 'breakpoint_calling')
    raw_data_directory = os.path.join(varcalls_dir, 'raw')
    breakpoints_filename = os.path.join(varcalls_dir, 'breakpoints.h5')
    ref_data_directory = config['ref_data_directory']

    pypeliner.workflow.Workflow(ctx={'docker_image': baseimage})

    workflow.setobj(
        obj=mgd.OutputChunks('tumour_cell_id'),
        value=tumour_bams.keys(),
    )

    if isinstance(normal_bams, dict):
        workflow.setobj(
            obj=mgd.OutputChunks('normal_cell_id'),
            value=normal_bams.keys(),
        )
        workflow.set_filenames('normal_cells.bam', 'normal_cell_id', fnames=normal_bams)
        normal_bam = mgd.InputFile('normal_cells.bam', 'normal_cell_id', extensions=['.bai'])
    else:
        normal_bam = mgd.InputFile(normal_bams, extensions=['.bai'])

    if args['destruct']:
        if isinstance(normal_bams, dict):
            workflow.transform(
                name='merge_normal_cells_destruct',
                func='single_cell.utils.bamutils.bam_merge',
                args=(
                    normal_bam,
                    mgd.TempOutputFile("merged_tumour.bam"),
                ),
                kwargs={'docker_image': config['docker']['samtools']}
            )
            final_wgs_bam = mgd.TempInputFile("merged_tumour.bam")
        else:
            final_wgs_bam = mgd.InputFile(normal_bams)

        workflow.subworkflow(
            name='destruct',
            ctx={'docker_image': baseimage},
            func="biowrappers.components.breakpoint_calling.destruct.destruct_pipeline",
            args=(
                final_wgs_bam,
                mgd.InputFile('tumour.bam', 'tumour_cell_id', fnames=tumour_bams),
                config.get('destruct', {}),
                ref_data_directory,
                mgd.OutputFile(breakpoints_filename),
                raw_data_directory,
            ),
        )

    if args['lumpy']:
        varcalls_dir = os.path.join(
            args['out_dir'], 'results', 'breakpoint_calling')
        breakpoints_bed = os.path.join(varcalls_dir, 'lumpy_breakpoints.bed')
        breakpoints_h5 = os.path.join(varcalls_dir, 'lumpy_breakpoints.h5')

        workflow.subworkflow(
            name='lumpy',
            ctx={'docker_image': baseimage},
            func="single_cell.workflows.lumpy.create_lumpy_workflow",
            args=(
                config,
                mgd.InputFile('tumour.bam', 'tumour_cell_id', fnames=tumour_bams, extensions=['.bai']),
                normal_bam,
                mgd.OutputFile(breakpoints_bed),
                mgd.OutputFile(breakpoints_h5),
                tumour_id,
                normal_id,
            ),
        )

    info_file = os.path.join(args["out_dir"],'results','breakpoint_calling', "info.yaml")

    results = {
        'destruct_data': helpers.format_file_yaml(breakpoints_filename),
    }

    tumour_dataset = {k: helpers.format_file_yaml(v) for k,v in tumour_bams.iteritems()}
    if isinstance(normal_bams, dict):
        normal_dataset = {k: helpers.format_file_yaml(v) for k, v in normal_bams.iteritems()}
    else:
        normal_dataset = helpers.format_file_yaml(normal_bams)
    input_datasets = {'normal': normal_dataset,
                      'tumour': tumour_dataset}

    metadata = {
        'breakpoint_calling': {
            'ref_data': ref_data_directory,
            'version': single_cell.__version__,
            'results': results,
            'containers': config['docker'],
            'input_datasets': input_datasets,
            'output_datasets': None
        }
    }

    workflow.transform(
        name='generate_meta_yaml',
        ctx=dict(mem=config['memory']['med'],
                 docker_image=baseimage,
                 mem_retry_increment=2, ncpus=1),
        func="single_cell.utils.helpers.write_to_yaml",
        args=(
            mgd.OutputFile(info_file),
            metadata
        )
    )


    return workflow

