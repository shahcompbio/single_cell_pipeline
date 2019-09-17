import pypeliner
import pypeliner.managed as mgd


def process_cells_destruct(
        destruct_config, cell_bam_files,
        reads_1, reads_2, sample_1, sample_2, stats,
        tag=False
):
    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1, }

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
        func="single_cell.workflows.destruct_singlecell.tasks.merge_read_counts",
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


def destruct_preprocess_workflow(
        normal_bam_files, normal_stats,
        normal_reads_1, normal_reads_2,
        normal_sample_1, normal_sample_2,
        ref_data_directory, destruct_config,
        tag=False
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
            value=list(normal_bam_files.keys()),
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
            kwargs={'tag': tag}
        )

    return workflow


def create_destruct_workflow(
        normal_stats, normal_reads_1, normal_reads_2, normal_sample_1, normal_sample_2,
        tumour_stats, tumour_reads_1, tumour_reads_2, tumour_sample_1, tumour_sample_2,
        destruct_config, ref_data_directory, breakpoints_filename,
        breakpoints_library_filename, cell_counts_filename, raw_data_directory,
        normal_sample_id='normal', tumour_sample_id='tumour',
        tumour_library_id='tumour'
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

    workflow.subworkflow(
        name='destruct',
        func="destruct.workflow.create_destruct_fastq_workflow",
        ctx={'disk': 200},
        args=(
            {
                normal_sample_id: mgd.InputFile(normal_reads_1),
                tumour_sample_id: mgd.InputFile(tumour_reads_1),
            },
            {
                normal_sample_id: mgd.InputFile(normal_reads_2),
                tumour_sample_id: mgd.InputFile(tumour_reads_2),
            },
            {
                normal_sample_id: mgd.InputFile(normal_sample_1),
                tumour_sample_id: mgd.InputFile(tumour_sample_1),
            },
            {
                normal_sample_id: mgd.InputFile(normal_sample_2),
                tumour_sample_id: mgd.InputFile(tumour_sample_2),
            },
            {
                normal_sample_id: mgd.InputFile(normal_stats),
                tumour_sample_id: mgd.InputFile(tumour_stats),
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
        name='filter_annotate_breakpoints',
        ctx={'mem': 8},
        func="biowrappers.components.breakpoint_calling.destruct.tasks.filter_annotate_breakpoints",
        args=(
            pypeliner.managed.TempInputFile('breakpoint_table'),
            pypeliner.managed.TempInputFile('breakpoint_library_table'),
            [normal_sample_id],
            pypeliner.managed.TempOutputFile('breakpoints_filename.csv'),
            pypeliner.managed.TempOutputFile('breakpoints_library_filename.csv'),
        ),
    )

    workflow.transform(
        name='filter_breakpoint_reads',
        ctx={'mem': 8},
        func="single_cell.workflows.destruct_singlecell.tasks.filter_reads_file",
        args=(
            mgd.TempInputFile('breakpoint_read_table'),
            pypeliner.managed.TempInputFile('breakpoints_filename.csv'),
            mgd.TempOutputFile('breakpoint_read_table_filtered'),
        ),
    )

    workflow.transform(
        name='extract_cell_counts',
        ctx={'mem': 8},
        func="single_cell.workflows.destruct_singlecell.tasks.extract_cell_counts",
        args=(
            mgd.TempInputFile('breakpoint_read_table_filtered'),
            mgd.TempOutputFile('cell_counts_filename.csv'),
        ),
    )

    workflow.transform(
        name='prep_cell_counts',
        ctx={'mem': 8, 'ncpus': 1},
        func="single_cell.utils.csvutils.prep_csv_files",
        args=(
            mgd.TempInputFile('cell_counts_filename.csv'),
            mgd.TempOutputFile("cell_counts_prep.csv.gz", extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_cell_counts',
        ctx={'mem': 8, 'ncpus': 1},
        func="single_cell.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("cell_counts_prep.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(cell_counts_filename, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='prep_breakpoints',
        ctx={'mem': 8, 'ncpus': 1},
        func="single_cell.utils.csvutils.prep_csv_files",
        args=(
            pypeliner.managed.TempInputFile('breakpoints_filename.csv'),
            pypeliner.managed.TempOutputFile("breakpoints_filename_prep.csv.gz", extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_breakpoints',
        ctx={'mem': 8, 'ncpus': 1},
        func="single_cell.utils.csvutils.finalize_csv",
        args=(
            pypeliner.managed.TempInputFile("breakpoints_filename_prep.csv.gz", extensions=['.yaml']),
            pypeliner.managed.OutputFile(breakpoints_filename, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='prep_breakpoints_library',
        ctx={'mem': 8, 'ncpus': 1},
        func="single_cell.utils.csvutils.prep_csv_files",
        args=(
            pypeliner.managed.TempInputFile('breakpoints_library_filename.csv'),
            pypeliner.managed.TempOutputFile('breakpoints_library_prep.csv.gz', extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_breakpoints_library',
        ctx={'mem': 8, 'ncpus': 1},
        func="single_cell.utils.csvutils.finalize_csv",
        args=(
            pypeliner.managed.TempInputFile('breakpoints_library_prep.csv.gz', extensions=['.yaml']),
            pypeliner.managed.OutputFile(breakpoints_library_filename, extensions=['.yaml']),
        ),
    )

    return workflow


def destruct_multi_sample_workflow(
        normal_bam, tumour_bam_files, destruct_config, config,
        destruct_ref_data_dir, breakpoints_csv, breakpoints_library_csv,
        cell_counts_csv, raw_data_dir, normal_sample_id='normal',
):
    ctx = {'docker_image': config['docker']['destruct']}
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'library_id', 'cell_id'),
        value=list(tumour_bam_files.keys()),
    )

    keys = [(sample_id, library_id) for (sample_id, library_id, _) in list(tumour_bam_files.keys())]
    keys = sorted(set(keys))

    breakpoints_csv = dict([(key, breakpoints_csv(*key))
                            for key in keys])
    breakpoints_library_csv = dict([(key, breakpoints_library_csv(*key))
                                    for key in keys])
    cell_counts_csv = dict([(key, cell_counts_csv(*key))
                            for key in keys])

    workflow.set_filenames('tumour_cells.bam', 'sample_id', 'library_id', 'cell_id', fnames=tumour_bam_files)
    workflow.set_filenames('breakpoints.csv', 'sample_id', 'library_id', fnames=breakpoints_csv)
    workflow.set_filenames('breakpoints_library.csv', 'sample_id', 'library_id', fnames=breakpoints_library_csv)
    workflow.set_filenames('cell_counts.csv', 'sample_id', 'library_id', fnames=cell_counts_csv)

    workflow.subworkflow(
        name='normal_preprocess_destruct',
        func='single_cell.workflows.destruct_singlecell.destruct_preprocess_workflow',
        args=(
            normal_bam,
            mgd.TempOutputFile('normal_stats'),
            mgd.TempOutputFile('normal_reads_1.fastq.gz'),
            mgd.TempOutputFile('normal_reads_2.fastq.gz'),
            mgd.TempOutputFile('normal_sample_1.fastq.gz'),
            mgd.TempOutputFile('normal_sample_2.fastq.gz'),
            destruct_ref_data_dir,
            destruct_config,
        ),
    )

    workflow.subworkflow(
        name='tumour_preprocess_destruct',
        func='single_cell.workflows.destruct_singlecell.destruct_preprocess_workflow',
        axes=('sample_id', 'library_id'),
        args=(
            mgd.InputFile('tumour_cells.bam', 'sample_id', 'library_id', 'cell_id', extensions=['.bai']),
            mgd.TempOutputFile('tumour_stats', 'sample_id', 'library_id'),
            mgd.TempOutputFile('tumour_reads_1.fastq.gz', 'sample_id', 'library_id'),
            mgd.TempOutputFile('tumour_reads_2.fastq.gz', 'sample_id', 'library_id'),
            mgd.TempOutputFile('tumour_sample_1.fastq.gz', 'sample_id', 'library_id'),
            mgd.TempOutputFile('tumour_sample_2.fastq.gz', 'sample_id', 'library_id'),
            destruct_ref_data_dir,
            destruct_config,
        ),
        kwargs={'tag': True}
    )

    workflow.subworkflow(
        name='run_destruct',
        func='single_cell.workflows.destruct_singlecell.create_destruct_workflow',
        axes=('sample_id', 'library_id'),
        args=(
            mgd.TempInputFile('normal_stats'),
            mgd.TempInputFile('normal_reads_1.fastq.gz'),
            mgd.TempInputFile('normal_reads_2.fastq.gz'),
            mgd.TempInputFile('normal_sample_1.fastq.gz'),
            mgd.TempInputFile('normal_sample_2.fastq.gz'),
            mgd.TempInputFile('tumour_stats', 'sample_id', 'library_id'),
            mgd.TempInputFile('tumour_reads_1.fastq.gz', 'sample_id', 'library_id'),
            mgd.TempInputFile('tumour_reads_2.fastq.gz', 'sample_id', 'library_id'),
            mgd.TempInputFile('tumour_sample_1.fastq.gz', 'sample_id', 'library_id'),
            mgd.TempInputFile('tumour_sample_2.fastq.gz', 'sample_id', 'library_id'),
            destruct_config,
            destruct_ref_data_dir,
            mgd.OutputFile('breakpoints.csv', 'sample_id', 'library_id'),
            mgd.OutputFile('breakpoints_library.csv', 'sample_id', 'library_id'),
            mgd.OutputFile('cell_counts.csv', 'sample_id', 'library_id'),
            mgd.Template(raw_data_dir, 'sample_id', 'library_id'),
        ),
        kwargs={
            'tumour_sample_id': mgd.Instance('sample_id'),
            'tumour_library_id': mgd.Instance('library_id'),
            'normal_sample_id': normal_sample_id,
        },
    )

    return workflow
