import pypeliner
import pypeliner.managed as mgd
from single_cell.workflows.destruct_singlecell.dtypes import dtypes


def process_cells_destruct(
        destruct_config, config, cell_bam_files,
        reads_1, reads_2, sample_1, sample_2, stats,
        tag=False
):

    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1, }

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    if isinstance(cell_bam_files, str):
        workflow.commandline(
            name='bamdisc_normal',
            # func="single_cell.workflows.destruct_singlecell.tasks.destruct_bamdisc_and_numreads",
            ctx={'io': 1, 'mem': 8, 'disk': 200, 'docker_image': config['docker']['destruct']},
            args=(
                'destruct_bamdiscordantfastq',
                '-r',
                '-c', destruct_config['bam_max_soft_clipped'],
                '-f', destruct_config['bam_max_fragment_length'],
                '-b', mgd.InputFile(cell_bam_files),
                '-s', mgd.OutputFile(stats),
                '--fastq1', mgd.OutputFile(reads_1),
                '--fastq2', mgd.OutputFile(reads_2),
                '-t', mgd.TempSpace('bamdisc_cell_tempspace'),
                '-n', destruct_config['num_read_samples'],
                '--sample1', mgd.OutputFile(sample_1),
                '--sample2', mgd.OutputFile(sample_2),
            )
        )

        return workflow

    cells = list(cell_bam_files.keys())

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cells,
    )

    workflow.commandline(
        name='bamdisc_and_numreads_cell',
        # func="single_cell.workflows.destruct_singlecell.tasks.destruct_bamdisc_and_numreads",
        axes=('cell_id',),
        ctx={'io': 1, 'mem': 8, 'docker_image': config['docker']['destruct']},
        args=(
            'destruct_bamdiscordantfastq',
            '-r',
            '-c', destruct_config['bam_max_soft_clipped'],
            '-f', destruct_config['bam_max_fragment_length'],
            '-b', mgd.InputFile('bam', 'cell_id', fnames=cell_bam_files),
            '-s', mgd.TempOutputFile('cell_stats', 'cell_id'),
            '--fastq1', mgd.TempOutputFile('cell_reads_1.fastq.gz', 'cell_id'),
            '--fastq2', mgd.TempOutputFile('cell_reads_2.fastq.gz', 'cell_id'),
            '-t', mgd.TempSpace('bamdisc_cell_tempspace', 'cell_id'),
            '-n', destruct_config['num_read_samples'],
            '--sample1', mgd.TempOutputFile('cell_sample_1.fastq.gz', 'cell_id'),
            '--sample2', mgd.TempOutputFile('cell_sample_2.fastq.gz', 'cell_id'),
        ),
    )

    workflow.transform(
        name='merge_reads_r1',
        ctx={'io': 1, 'mem': 8, 'disk': 100},
        func="single_cell.workflows.destruct_singlecell.tasks.merge_fastqs",
        args=(
            mgd.TempInputFile('cell_reads_1.fastq.gz', 'cell_id'),
            mgd.OutputFile(reads_1),
        ),
        kwargs={'tag': tag}
    )

    workflow.transform(
        name='merge_reads_r2',
        ctx={'io': 1, 'mem': 8, 'disk': 100},
        func="single_cell.workflows.destruct_singlecell.tasks.merge_fastqs",
        args=(
            mgd.TempInputFile('cell_reads_2.fastq.gz', 'cell_id'),
            mgd.OutputFile(reads_2),
        ),
        kwargs={'tag': tag}
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
        ref_data_directory, destruct_config, config,
        tag=False
):
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name="get_destruct_config",
        func="destruct.defaultconfig.get_config",
        ctx={'docker_image': config['docker']['destruct'], 'disk': 200},
        ret=mgd.TempOutputObj("destruct_config"),
        args=(
            ref_data_directory,
            destruct_config
        )
    )

    if isinstance(normal_bam_files, str):
        workflow.subworkflow(
            name='process_individual_cells',
            func=process_cells_destruct,
            args=(
                mgd.TempInputObj("destruct_config"),
                config,
                mgd.InputFile(normal_bam_files),
                mgd.OutputFile(normal_reads_1),
                mgd.OutputFile(normal_reads_2),
                mgd.OutputFile(normal_sample_1),
                mgd.OutputFile(normal_sample_2),
                mgd.OutputFile(normal_stats),
            ),
            kwargs={'tag': tag}
        )
    else:
        workflow.setobj(
            obj=mgd.OutputChunks('normal_cell_id'),
            value=list(normal_bam_files.keys()),
        )

        workflow.subworkflow(
            name='process_individual_cells',
            func=process_cells_destruct,
            args=(
                mgd.TempInputObj("destruct_config"),
                config,
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


def destruct_workflow(
        normal_stats, normal_reads_1, normal_reads_2, normal_sample_1, normal_sample_2,
        tumour_stats, tumour_reads_1, tumour_reads_2, tumour_sample_1, tumour_sample_2,
        destruct_config, config, ref_data_directory, breakpoints_filename,
        breakpoints_library_filename, cell_counts_filename, raw_data_directory,
        normal_sample_id='normal', tumour_sample_id='tumour',
        tumour_library_id='tumour'
):
    tumour_sample_id = '_'.join([tumour_sample_id, tumour_library_id])
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name="get_destruct_config",
        func="destruct.defaultconfig.get_config",
        ctx={'docker_image': config['docker']['destruct']},
        ret=mgd.TempOutputObj("destruct_config"),
        args=(
            ref_data_directory,
            destruct_config
        )
    )

    workflow.subworkflow(
        name='destruct',
        func="destruct.workflow.create_destruct_fastq_workflow",
        ctx={'docker_image': config['docker']['destruct'], 'disk': 200},
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
        ctx={'docker_image': config['docker']['destruct'], 'mem': 8},
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
        func="single_cell.utils.csvutils.rewrite_csv_file",
        args=(
            mgd.TempInputFile('cell_counts_filename.csv'),
            mgd.OutputFile(cell_counts_filename, extensions=['.yaml']),
        ),
        kwargs={'dtypes': dtypes()['cell_counts']}
    )

    workflow.transform(
        name='prep_breakpoints',
        ctx={'mem': 8, 'ncpus': 1},
        func="single_cell.utils.csvutils.rewrite_csv_file",
        args=(
            pypeliner.managed.TempInputFile('breakpoints_filename.csv'),
            pypeliner.managed.OutputFile(breakpoints_filename, extensions=['.yaml']),
        ),
        kwargs={'dtypes': dtypes()['breakpoints']}
    )


    workflow.transform(
        name='prep_breakpoints_library',
        ctx={'mem': 8, 'ncpus': 1},
        func="single_cell.utils.csvutils.rewrite_csv_file",
        args=(
            pypeliner.managed.TempInputFile('breakpoints_library_filename.csv'),
            pypeliner.managed.OutputFile(breakpoints_library_filename, extensions=['.yaml']),
        ),
        kwargs={'dtypes': dtypes()['library']}
    )

    return workflow


def create_destruct_workflow(
        normal_bam, tumour_bam_files, destruct_config, config,
        destruct_ref_data_dir, breakpoints_csv, breakpoints_library_csv,
        cell_counts_csv, normal_sample_id='normal',
):
    ctx = {'docker_image': config['docker']['single_cell_pipeline']}
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=list(tumour_bam_files.keys()),
    )

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
            config
        ),
    )

    workflow.subworkflow(
        name='tumour_preprocess_destruct',
        func='single_cell.workflows.destruct_singlecell.destruct_preprocess_workflow',
        args=(
            mgd.InputFile('tumour_cells.bam', 'cell_id', extensions=['.bai'], fnames=tumour_bam_files),
            mgd.TempOutputFile('tumour_stats'),
            mgd.TempOutputFile('tumour_reads_1.fastq.gz'),
            mgd.TempOutputFile('tumour_reads_2.fastq.gz'),
            mgd.TempOutputFile('tumour_sample_1.fastq.gz'),
            mgd.TempOutputFile('tumour_sample_2.fastq.gz'),
            destruct_ref_data_dir,
            destruct_config,
            config
        ),
        kwargs={'tag': True}
    )

    workflow.subworkflow(
        name='run_destruct',
        func='single_cell.workflows.destruct_singlecell.destruct_workflow',
        args=(
            mgd.TempInputFile('normal_stats'),
            mgd.TempInputFile('normal_reads_1.fastq.gz'),
            mgd.TempInputFile('normal_reads_2.fastq.gz'),
            mgd.TempInputFile('normal_sample_1.fastq.gz'),
            mgd.TempInputFile('normal_sample_2.fastq.gz'),
            mgd.TempInputFile('tumour_stats'),
            mgd.TempInputFile('tumour_reads_1.fastq.gz'),
            mgd.TempInputFile('tumour_reads_2.fastq.gz'),
            mgd.TempInputFile('tumour_sample_1.fastq.gz'),
            mgd.TempInputFile('tumour_sample_2.fastq.gz'),
            destruct_config,
            config,
            destruct_ref_data_dir,
            mgd.OutputFile(breakpoints_csv),
            mgd.OutputFile(breakpoints_library_csv),
            mgd.OutputFile(cell_counts_csv),
            mgd.TempSpace("raw_data_dir"),
        ),
    )

    return workflow
