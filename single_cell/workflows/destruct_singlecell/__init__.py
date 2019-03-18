import pypeliner
import pypeliner.managed as mgd


def process_cells_destruct(
        destruct_config, cell_bam_files,
        reads_1, reads_2, sample_1, sample_2, stats,
        tag=False
):

    ctx = {'mem_retry_increment': 2, 'ncpus': 1,}

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cell_bam_files.keys(),
    )

    workflow.transform(
        name='bamdisc_cell',
        func="single_cell.workflows.destruct_singlecell.tasks.destruct_bamdisc",
        axes=('cell_id',),
        ctx={'io': 1, 'mem': 8},
        args=(
            destruct_config,
            mgd.InputFile('bam', 'cell_id', fnames=cell_bam_files),
            mgd.TempOutputFile('cell_stats', 'cell_id'),
            mgd.TempOutputFile('cell_reads_1.fastq.gz', 'cell_id'),
            mgd.TempOutputFile('cell_reads_2.fastq.gz', 'cell_id'),
            mgd.TempOutputFile('cell_sample_1.fastq.gz', 'cell_id'),
            mgd.TempOutputFile('cell_sample_2.fastq.gz', 'cell_id'),
            mgd.TempSpace('bamdisc_cell_tempspace', 'cell_id'),
        )
    )

    workflow.transform(
        name='merge_reads',
        ctx={'io': 1, 'mem': 8},
        func="single_cell.workflows.destruct_singlecell.tasks.merge_cell_fastqs",
        args=(
            mgd.TempInputFile('cell_reads_1.fastq.gz', 'cell_id'),
            mgd.TempInputFile('cell_reads_2.fastq.gz', 'cell_id'),
            mgd.OutputFile(reads_1),
            mgd.OutputFile(reads_2),
        ),
        kwargs={'tag': tag}
    )

    workflow.transform(
        name='merge_sample',
        ctx={'io': 1, 'mem': 8},
        func="single_cell.workflows.destruct_singlecell.tasks.resample_fastqs",
        args=(
            mgd.TempInputFile('cell_sample_1.fastq.gz', 'cell_id'),
            mgd.TempInputFile('cell_sample_2.fastq.gz', 'cell_id'),
            mgd.OutputFile(sample_1),
            mgd.OutputFile(sample_2),
            destruct_config['num_read_samples'],
        ),
    )

    workflow.transform(
        name='merge_stats',
        ctx={'io': 1, 'mem': 8},
        func="single_cell.workflows.destruct_singlecell.tasks.merge_stats",
        args=(
            mgd.TempInputFile('cell_stats', 'cell_id'),
            mgd.OutputFile(stats),
        ),
    )

    return workflow


def create_destruct_workflow(
    normal_bam_files,
    tumour_bam_files,
    destruct_config,
    ref_data_directory,
    breakpoints_filename,
    breakpoints_library_filename,
    cell_counts_filename,
    raw_data_directory,
    normal_sample_id='normal',
    tumour_sample_id='tumour',
):
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name="get_destruct_config",
        func="destruct.defaultconfig.get_config",
        ret=mgd.TempOutputObj("destruct_config"),
        args=(
            ref_data_directory,
            destruct_config
        )
    )

    if isinstance(normal_bam_files, str):
        workflow.transform(
            name='bamdisc_normal',
            func="single_cell.workflows.destruct_singlecell.tasks.destruct_bamdisc",
            ctx={'io': 1, 'mem': 8},
            args=(
                mgd.TempInputObj("destruct_config"),
                mgd.InputFile(normal_bam_files),
                mgd.TempOutputFile('normal_stats'),
                mgd.TempOutputFile('normal_reads_1.fastq.gz'),
                mgd.TempOutputFile('normal_reads_2.fastq.gz'),
                mgd.TempOutputFile('normal_sample_1.fastq.gz'),
                mgd.TempOutputFile('normal_sample_2.fastq.gz'),
                mgd.TempSpace('bamdisc_normal_tempspace'),
            )
        )
    else:
        workflow.setobj(
            obj=mgd.OutputChunks('normal_cell_id'),
            value=normal_bam_files.keys(),
        )

        workflow.subworkflow(
            name='process_normal_cells',
            func=process_cells_destruct,
            args=(
                mgd.TempInputObj("destruct_config"),
                mgd.InputFile('bam', 'normal_cell_id', fnames=normal_bam_files),
                mgd.TempOutputFile('normal_reads_1.fastq.gz'),
                mgd.TempOutputFile('normal_reads_2.fastq.gz'),
                mgd.TempOutputFile('normal_sample_1.fastq.gz'),
                mgd.TempOutputFile('normal_sample_2.fastq.gz'),
                mgd.TempOutputFile('normal_stats'),
            ),
            kwargs={'tag': False}
        )

    workflow.setobj(
        obj=mgd.OutputChunks('tumour_cell_id'),
        value=tumour_bam_files.keys(),
    )

    workflow.subworkflow(
        name='process_tumour_cells',
        func=process_cells_destruct,
        args=(
            mgd.TempInputObj("destruct_config"),
            mgd.InputFile('bam', 'tumour_cell_id', fnames=tumour_bam_files),
            mgd.TempOutputFile('tumour_reads_1.fastq.gz'),
            mgd.TempOutputFile('tumour_reads_2.fastq.gz'),
            mgd.TempOutputFile('tumour_sample_1.fastq.gz'),
            mgd.TempOutputFile('tumour_sample_2.fastq.gz'),
            mgd.TempOutputFile('tumour_stats'),
        ),
        kwargs={'tag': True}
    )

    workflow.subworkflow(
        name='destruct',
        func="destruct.workflow.create_destruct_fastq_workflow",
        args=(
            {
                normal_sample_id: mgd.TempInputFile('normal_reads_1.fastq.gz'),
                tumour_sample_id: mgd.TempInputFile('tumour_reads_1.fastq.gz'),
            },
            {
                normal_sample_id: mgd.TempInputFile('normal_reads_2.fastq.gz'),
                tumour_sample_id: mgd.TempInputFile('tumour_reads_2.fastq.gz'),
            },
            {
                normal_sample_id: mgd.TempInputFile('normal_sample_1.fastq.gz'),
                tumour_sample_id: mgd.TempInputFile('tumour_sample_1.fastq.gz'),
            },
            {
                normal_sample_id: mgd.TempInputFile('normal_sample_2.fastq.gz'),
                tumour_sample_id: mgd.TempInputFile('tumour_sample_2.fastq.gz'),
            },
            {
                normal_sample_id: mgd.TempInputFile('normal_stats'),
                tumour_sample_id: mgd.TempInputFile('tumour_stats'),
            },
            mgd.TempOutputFile('breakpoint_table'),
            mgd.TempOutputFile('breakpoint_library_table'),
            mgd.TempOutputFile('breakpoint_read_table'),
            mgd.TempInputObj("destruct_config"),
            ref_data_directory,
        ),
        kwargs={
            'raw_data_dir': raw_data_directory,
        },
    )

    workflow.transform(
        name='extract_cell_counts',
        ctx={'mem': 8},
        func="single_cell.workflows.destruct_singlecell.tasks.extract_cell_counts",
        args=(
            mgd.TempInputFile('breakpoint_read_table'),
            mgd.OutputFile(cell_counts_filename),
        ),
    )

    workflow.transform(
        name='filter_annotate_breakpoints',
        ctx={'mem': 8},
        func="biowrappers.components.breakpoint_calling.destruct.tasks.filter_annotate_breakpoints",
        args=(
            pypeliner.managed.TempInputFile('breakpoint_table'),
            pypeliner.managed.TempInputFile('breakpoint_library_table'),
            [normal_sample_id],
            pypeliner.managed.OutputFile(breakpoints_filename),
            pypeliner.managed.OutputFile(breakpoints_library_filename),
        ),
    )

    return workflow

