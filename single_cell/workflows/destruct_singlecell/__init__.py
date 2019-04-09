import pypeliner
import pypeliner.managed as mgd


def process_cells_destruct(
        destruct_config, cell_bam_files,
        reads_1, reads_2, sample_1, sample_2, stats,
        tag=False
):

    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1,}

    cells = list(cell_bam_files.keys())

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cells,
    )

    workflow.transform(
        name='bamdisc_and_numreads_cell',
        func="single_cell.workflows.destruct_singlecell.tasks.destruct_bamdisc_and_numreads",
        axes=('cell_id',),
        ctx={'io': 1, 'mem': 8},
        ret=mgd.TempOutputObj("numreads", "cell_id"),
        args=(
            destruct_config,
            mgd.InputFile('bam', 'cell_id', fnames=cell_bam_files),
            mgd.TempOutputFile('cell_stats', 'cell_id'),
            mgd.TempOutputFile('cell_reads_1.fastq.gz', 'cell_id'),
            mgd.TempOutputFile('cell_reads_2.fastq.gz', 'cell_id'),
            mgd.TempOutputFile('cell_sample_1.fastq.gz', 'cell_id'),
            mgd.TempOutputFile('cell_sample_2.fastq.gz', 'cell_id'),
            mgd.TempSpace('bamdisc_cell_tempspace', 'cell_id'),
        ),
    )

    workflow.transform(
        name='merge_read_counts',
        ret=mgd.TempOutputObj("readcounts"),
        func="single_cell.workflows.destruct_singlecell.tasks.get_max_read_count",
        ctx={'io': 1, 'mem': 8},
        args=(
            mgd.TempInputObj('numreads', 'cell_id'),
        )
    )

    workflow.transform(
        name='reindex_reads',
        func="single_cell.workflows.destruct_singlecell.tasks.re_index_reads_both",
        ctx={'io': 1, 'mem': 8},
        axes=('cell_id',),
        args=(
            mgd.TempInputFile('cell_reads_1.fastq.gz', 'cell_id'),
            mgd.TempOutputFile('cell_reads_1_reindex.fastq.gz', 'cell_id'),
            mgd.TempInputFile('cell_reads_2.fastq.gz', 'cell_id'),
            mgd.TempOutputFile('cell_reads_2_reindex.fastq.gz', 'cell_id'),
            mgd.InputInstance('cell_id'),
            cells,
            mgd.TempInputObj('readcounts'),
        ),
        kwargs={'tag': tag}
    )

    workflow.transform(
        name='merge_reads_r1',
        ctx={'io': 1, 'mem': 8, 'disk': 100},
        func="single_cell.workflows.destruct_singlecell.tasks.merge_cell_fastqs",
        args=(
            mgd.TempInputFile('cell_reads_1_reindex.fastq.gz', 'cell_id'),
            mgd.OutputFile(reads_1),
        ),
    )

    workflow.transform(
        name='merge_reads_r2',
        ctx={'io': 1, 'mem': 8, 'disk': 100},
        func="single_cell.workflows.destruct_singlecell.tasks.merge_cell_fastqs",
        args=(
            mgd.TempInputFile('cell_reads_2_reindex.fastq.gz', 'cell_id'),
            mgd.OutputFile(reads_2),
        ),
    )

    workflow.transform(
        name='merge_sample',
        ctx={'io': 1, 'mem': 8, 'disk': 100},
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


def destruct_normal_preprocess_workflow(
        normal_bam_files, normal_stats,
        normal_reads_1, normal_reads_2,
        normal_sample_1, normal_sample_2,
        ref_data_directory, destruct_config
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
            func="single_cell.workflows.destruct_singlecell.tasks.destruct_bamdisc_and_numreads",
            ctx={'io': 1, 'mem': 8, 'disk': 200},
            args=(
                mgd.TempInputObj("destruct_config"),
                mgd.InputFile(normal_bam_files),
                mgd.OutputFile(normal_stats),
                mgd.OutputFile(normal_reads_1),
                mgd.OutputFile(normal_reads_2),
                mgd.OutputFile(normal_sample_1),
                mgd.OutputFile(normal_sample_2),
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
                mgd.OutputFile(normal_reads_1),
                mgd.OutputFile(normal_reads_2),
                mgd.OutputFile(normal_sample_1),
                mgd.OutputFile(normal_sample_2),
                mgd.OutputFile(normal_stats),
            ),
            kwargs={'tag': False}
        )

    return workflow


def create_destruct_workflow(
    tumour_bam_files,
    normal_stats,
    normal_reads_1,
    normal_reads_2,
    normal_sample_1,
    normal_sample_2,
    destruct_config,
    ref_data_directory,
    breakpoints_filename,
    breakpoints_library_filename,
    cell_counts_filename,
    raw_data_directory,
    normal_sample_id='normal',
    tumour_sample_id='tumour',
    tumour_library_id='tumour',
):
    tumour_sample_id = '_'.join([tumour_sample_id, tumour_library_id])
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

    workflow.setobj(
        obj=mgd.OutputChunks('tumour_cell_id'),
        value=list(tumour_bam_files.keys()),
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
                normal_sample_id: mgd.InputFile(normal_reads_1),
                tumour_sample_id: mgd.TempInputFile('tumour_reads_1.fastq.gz'),
            },
            {
                normal_sample_id: mgd.InputFile(normal_reads_2),
                tumour_sample_id: mgd.TempInputFile('tumour_reads_2.fastq.gz'),
            },
            {
                normal_sample_id: mgd.InputFile(normal_sample_1),
                tumour_sample_id: mgd.TempInputFile('tumour_sample_1.fastq.gz'),
            },
            {
                normal_sample_id: mgd.InputFile(normal_sample_2),
                tumour_sample_id: mgd.TempInputFile('tumour_sample_2.fastq.gz'),
            },
            {
                normal_sample_id: mgd.InputFile(normal_stats),
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

