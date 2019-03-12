import os
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers
import single_cell


def create_destruct_workflow(
    normal_bam_file,
    tumour_bam_files,
    config,
    ref_data_directory,
    breakpoints_filename,
    raw_data_directory,
    normal_sample_id='normal',
):
    bam_files = tumour_bam_files
    bam_files[normal_sample_id] = normal_bam_file

    baseimage = config['docker']['single_cell_pipeline']

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=bam_files.keys(),
    )

    workflow.commandline(
        name='bamdisc',
        axes=('cell_id',),
        ctx={'io': 1, 'mem': 8, 'docker_image': baseimage},
        args=(
            'destruct_bamdiscordantfastq',
            '-r',
            '-c', config['breakpoint_calling']['bam_max_soft_clipped'],
            '-f', config['breakpoint_calling']['bam_max_fragment_length'],
            '-b', mgd.InputFile('bam', 'cell_id', fnames=bam_files),
            '-s', mgd.TempOutputFile('stats', 'cell_id'),
            '--fastq1', mgd.TempOutputFile('reads_1', 'cell_id'),
            '--fastq2', mgd.TempOutputFile('reads_2', 'cell_id'),
            '-t', mgd.TempSpace('bamdisc_tempspace', 'cell_id'),
            '-n', config['breakpoint_calling']['num_read_samples'],
            '--sample1', mgd.TempOutputFile('sample_1', 'cell_id'),
            '--sample2', mgd.TempOutputFile('sample_2', 'cell_id'),
        ),
    )

    workflow.transform(
        name='merge_reads',
        ctx={'io': 1, 'mem': 8},
        func="single_cell.workflows.destruct.tasks.merge_fastqs",
        args=(
            mgd.TempInputFile('reads_1', 'cell_id'),
            mgd.TempInputFile('reads_2', 'cell_id'),
            mgd.TempOutputFile('reads_1'),
            mgd.TempOutputFile('reads_2'),
        ),
    )

    workflow.transform(
        name='merge_sample',
        ctx={'io': 1, 'mem': 8},
        func="single_cell.workflows.destruct.tasks.resample_fastqs",
        args=(
            mgd.TempInputFile('sample_1', 'cell_id'),
            mgd.TempInputFile('sample_2', 'cell_id'),
            mgd.TempOutputFile('sample_1'),
            mgd.TempOutputFile('sample_2'),
            config['breakpoint_calling']['num_read_samples'],
        ),
    )

    workflow.transform(
        name='merge_stats',
        ctx={'io': 1, 'mem': 8},
        func="single_cell.workflows.destruct.tasks.merge_stats",
        args=(
            mgd.TempInputFile('stats', 'cell_id'),
            mgd.TempOutputFile('stats'),
        ),
    )

    workflow.subworkflow(
        name='destruct',
        ctx={'docker_image': baseimage},
        func="destruct.workflows.create_destruct_fastq_workflow",
        args=(
            mgd.TempInputFile('reads_1'),
            mgd.TempInputFile('reads_2'),
            mgd.TempInputFile('sample_1'),
            mgd.TempInputFile('sample_2'),
            mgd.TempInputFile('stats'),
            mgd.TempOutputFile('breakpoint_table'),
            mgd.TempOutputFile('breakpoint_library_table'),
            mgd.TempOutputFile('breakpoint_read_table'),
            config['breakpoint_calling']['destruct'],
            ref_data_directory,
        ),
        kwargs={
            'raw_data_dir': raw_data_dir,
        },
    )




    workflow.transform(
        name='filter_annotate_breakpoints',
        ctx={'mem': 8},
        func=tasks.filter_annotate_breakpoints,
        args=(
            pypeliner.managed.InputFile(breakpoint_file),
            pypeliner.managed.InputFile(breakpoint_library_file),
            [normal_sample_id],
            pypeliner.managed.OutputFile(somatic_breakpoint_file),
            pypeliner.managed.OutputFile(somatic_breakpoint_library_file),
        ),
    )



    return workflow

