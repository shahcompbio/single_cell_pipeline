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
    tumour_sample_id='tumour',
):
    baseimage = config['docker']['single_cell_pipeline']

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=tumour_bam_files.keys(),
    )

    workflow.commandline(
        name='bamdisc_normal',
        ctx={'io': 1, 'mem': 8, 'docker_image': baseimage},
        args=(
            'destruct_bamdiscordantfastq',
            '-r',
            '-c', config['breakpoint_calling']['bam_max_soft_clipped'],
            '-f', config['breakpoint_calling']['bam_max_fragment_length'],
            '-b', mgd.InputFile(normal_bam_file),
            '-s', mgd.TempOutputFile('normal_stats'),
            '--fastq1', mgd.TempOutputFile('normal_reads_1'),
            '--fastq2', mgd.TempOutputFile('normal_reads_2'),
            '-t', mgd.TempSpace('bamdisc_normal_tempspace'),
            '-n', config['breakpoint_calling']['num_read_samples'],
            '--sample1', mgd.TempOutputFile('normal_sample_1'),
            '--sample2', mgd.TempOutputFile('normal_sample_2'),
        ),
    )

    workflow.commandline(
        name='bamdisc_tumour',
        axes=('cell_id',),
        ctx={'io': 1, 'mem': 8, 'docker_image': baseimage},
        args=(
            'destruct_bamdiscordantfastq',
            '-r',
            '-c', config['breakpoint_calling']['bam_max_soft_clipped'],
            '-f', config['breakpoint_calling']['bam_max_fragment_length'],
            '-b', mgd.InputFile('bam', 'cell_id', fnames=bam_files),
            '-s', mgd.TempOutputFile('tumour_stats', 'cell_id'),
            '--fastq1', mgd.TempOutputFile('tumour_reads_1', 'cell_id'),
            '--fastq2', mgd.TempOutputFile('tumour_reads_2', 'cell_id'),
            '-t', mgd.TempSpace('bamdisc_tumour_tempspace', 'cell_id'),
            '-n', config['breakpoint_calling']['num_read_samples'],
            '--sample1', mgd.TempOutputFile('tumour_sample_1', 'cell_id'),
            '--sample2', mgd.TempOutputFile('tumour_sample_2', 'cell_id'),
        ),
    )

    workflow.transform(
        name='merge_reads',
        ctx={'io': 1, 'mem': 8},
        func="single_cell.workflows.destruct.tasks.merge_fastqs",
        args=(
            mgd.TempInputFile('tumour_reads_1', 'cell_id'),
            mgd.TempInputFile('tumour_reads_2', 'cell_id'),
            mgd.TempOutputFile('tumour_reads_1'),
            mgd.TempOutputFile('tumour_reads_2'),
        ),
    )

    workflow.transform(
        name='merge_sample',
        ctx={'io': 1, 'mem': 8},
        func="single_cell.workflows.destruct.tasks.resample_fastqs",
        args=(
            mgd.TempInputFile('tumour_sample_1', 'cell_id'),
            mgd.TempInputFile('tumour_sample_2', 'cell_id'),
            mgd.TempOutputFile('tumour_sample_1'),
            mgd.TempOutputFile('tumour_sample_2'),
            config['breakpoint_calling']['num_read_samples'],
        ),
    )

    workflow.transform(
        name='merge_stats',
        ctx={'io': 1, 'mem': 8},
        func="single_cell.workflows.destruct.tasks.merge_stats",
        args=(
            mgd.TempInputFile('tumour_stats', 'cell_id'),
            mgd.TempOutputFile('tumour_stats'),
        ),
    )

    workflow.subworkflow(
        name='destruct',
        ctx={'docker_image': baseimage},
        func="destruct.workflows.create_destruct_fastq_workflow",
        args=(
            {
                normal_sample_id: mgd.TempInputFile('normal_reads_1'),
                tumour_sample_id: mgd.TempInputFile('tumour_reads_1'),
            },
            {
                normal_sample_id: mgd.TempInputFile('normal_reads_2'),
                tumour_sample_id: mgd.TempInputFile('tumour_reads_2'),
            },
            {
                normal_sample_id: mgd.TempInputFile('normal_sample_1'),
                tumour_sample_id: mgd.TempInputFile('tumour_sample_1'),
            },
            {
                normal_sample_id: mgd.TempInputFile('normal_sample_2'),
                tumour_sample_id: mgd.TempInputFile('tumour_sample_2'),
            },
            {
                normal_sample_id: mgd.TempInputFile('normal_stats'),
                tumour_sample_id: mgd.TempInputFile('tumour_stats'),
            },
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
        name='extract_cell_counts',
        ctx={'mem': 8},
        func=
        args=(
            mgd.TempInputFile('breakpoint_read_table'),
            mgd.TempOutputFile('cell_counts'),
        ),
    )

    workflow.transform(
        name='filter_annotate_breakpoints',
        ctx={'mem': 8},
        func=tasks.filter_annotate_breakpoints,
        args=(
            pypeliner.managed.TempInputFile('breakpoint_table'),
            pypeliner.managed.TempInputFile('breakpoint_library_table'),
            pypeliner.managed.TempInputFile('cell_counts'),
            [normal_sample_id],
            pypeliner.managed.OutputFile(breakpoints_filename),
        ),
    )

    return workflow

