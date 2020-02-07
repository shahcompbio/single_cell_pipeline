import pypeliner
import pypeliner.managed as mgd
from single_cell.workflows.lumpy.dtypes import dtypes


def lumpy_preprocess_workflow(
        bam_files, config, discordants, split_reads,
        histogram, mean_stdev
):
    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1,
           'docker_image': config['docker']['single_cell_pipeline']
           }

    lumpydocker = {'docker_image': config['docker']['lumpy']}

    histogram_settings = dict(
        N=10000, skip=0, min_elements=100, mads=10, X=4, read_length=101
    )
    histogram_settings.update(lumpydocker)

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    if isinstance(bam_files, dict):
        workflow.setobj(
            obj=mgd.OutputChunks('cell_id'),
            value=list(bam_files.keys()),
        )
        workflow.set_filenames('cells.bam', 'cell_id', fnames=bam_files)
        workflow.subworkflow(
            name='process_cells',
            func='single_cell.workflows.lumpy.lumpy_preprocess_cells',
            args=(
                config,
                mgd.InputFile('cells.bam', 'cell_id', fnames=bam_files, extensions=['.bai']),
                mgd.OutputFile(discordants),
                mgd.OutputFile(split_reads),
                mgd.OutputFile(histogram),
                mgd.OutputFile(mean_stdev)
            ),
        )
    else:
        workflow.transform(
            name='process_bulk',
            ctx={'mem': 8, 'ncpus': 1, 'disk': 200},
            func='single_cell.workflows.lumpy.tasks.process_bam',
            args=(
                mgd.InputFile(bam_files, extensions=['.bai']),
                mgd.OutputFile(discordants),
                mgd.OutputFile(split_reads),
                mgd.TempOutputFile('hist_normal.csv'),
                mgd.TempSpace("lumpy_normal_processing"),
            ),
            kwargs=histogram_settings,
        )
        workflow.transform(
            name='format_histo_bulk',
            ctx={'mem': 8, 'ncpus': 1},
            func='single_cell.workflows.lumpy.merge_histograms.merge_histograms',
            args=(
                mgd.TempInputFile('hist_normal.csv'),
                mgd.OutputFile(histogram),
                mgd.OutputFile(mean_stdev)
            ),
        )

    return workflow


def lumpy_preprocess_cells(
        config, bam_files, merged_discordants, merged_splitters, hist_csv, mean_stdev_obj
):
    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1,
           'docker_image': config['docker']['single_cell_pipeline']
           }

    lumpydocker = {'docker_image': config['docker']['lumpy']}

    histogram_settings = dict(
        N=10000, skip=0, min_elements=100, mads=10, X=4, read_length=101
    )
    histogram_settings.update(lumpydocker)

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=list(bam_files.keys()),
    )

    workflow.transform(
        name='process_tumour_cells',
        axes=('cell_id',),
        ctx={'mem': 8, 'ncpus': 1},
        func='single_cell.workflows.lumpy.tasks.process_bam',
        args=(
            mgd.InputFile('tumour_bam', 'cell_id', fnames=bam_files, extensions=['.bai']),
            mgd.TempOutputFile('tumour.discordants.sorted.bam', 'cell_id'),
            mgd.TempOutputFile('tumour.splitters.sorted.bam', 'cell_id'),
            mgd.TempOutputFile('hist.csv', 'cell_id'),
            mgd.TempSpace("lumpy_tumour_processing", "cell_id"),
        ),
        kwargs=dict(tag=mgd.InputInstance('cell_id'), **histogram_settings),
    )

    workflow.transform(
        name='merge_disc',
        ctx={'mem': 8, 'ncpus': 1},
        func='single_cell.workflows.lumpy.tasks.merge_bams',
        args=(
            mgd.TempInputFile('tumour.discordants.sorted.bam', 'cell_id'),
            mgd.OutputFile(merged_discordants),
            mgd.TempSpace("merge_disc_temp")
        ),
        kwargs=lumpydocker,
    )

    workflow.transform(
        name='merge_split',
        ctx={'mem': 8, 'ncpus': 1},
        func='single_cell.workflows.lumpy.tasks.merge_bams',
        args=(
            mgd.TempInputFile('tumour.splitters.sorted.bam', 'cell_id'),
            mgd.OutputFile(merged_splitters),
            mgd.TempSpace("merge_split_temp")
        ),
        kwargs=lumpydocker,
    )

    workflow.transform(
        name='merge_histo',
        ctx={'mem': 8, 'ncpus': 1},
        func='single_cell.workflows.lumpy.merge_histograms.merge_histograms',
        args=(
            mgd.TempInputFile('hist.csv', 'cell_id'),
            mgd.OutputFile(hist_csv),
            mgd.OutputFile(mean_stdev_obj)
        ),
    )

    return workflow


def lumpy_calling_workflow(
        config,
        normal_disc_reads, normal_split_reads, normal_histogram, normal_mean_stdev,
        tumour_disc_reads, tumour_split_reads, tumour_histogram, tumour_mean_stdev,
        lumpy_bed, lumpy_calls, lumpy_evidence, tumour_id='tumour', normal_id='normal',
        sample_id=None, library_id=None
):
    if sample_id and not tumour_id:
        tumour_id = sample_id + '_' + library_id
    if sample_id and not normal_id:
        normal_id = sample_id + 'N'

    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1,
           'docker_image': config['docker']['single_cell_pipeline']
           }

    lumpydocker = {'docker_image': config['docker']['lumpy']}

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.transform(
        name='run_lumpy',
        ctx={'mem': 8, 'ncpus': 1},
        func='single_cell.workflows.lumpy.tasks.run_lumpy',
        args=(
            mgd.InputFile(tumour_disc_reads),
            mgd.InputFile(tumour_split_reads),
            mgd.InputFile(tumour_histogram),
            mgd.InputFile(tumour_mean_stdev),
            tumour_id,
            mgd.InputFile(normal_disc_reads),
            mgd.InputFile(normal_split_reads),
            mgd.InputFile(normal_histogram),
            mgd.InputFile(normal_mean_stdev),
            normal_id,
            mgd.OutputFile(lumpy_bed),
            mgd.TempSpace("lumpy_temp"),
        ),
        kwargs=lumpydocker,
    )

    workflow.transform(
        name='parse_lumpy',
        ctx={'mem': 8, 'ncpus': 1},
        func='single_cell.workflows.lumpy.parse_lumpy_to_csv.parse_lumpy',
        args=(
            mgd.InputFile(lumpy_bed),
            mgd.TempOutputFile('lumpy_calls.csv.gz'),
            mgd.TempOutputFile('lumpy_evidence.csv.gz'),
        ),
    )

    workflow.transform(
        name='prep_lumpy_calls',
        ctx={'mem': 8, 'ncpus': 1},
        func="single_cell.utils.csvutils.rewrite_csv_file",
        args=(
            mgd.TempInputFile('lumpy_calls.csv.gz'),
            mgd.OutputFile(lumpy_calls, extensions=['.yaml']),
        ),
        kwargs= {
            'dtypes': dtypes()['breakpoint']
        }
    )

    workflow.transform(
        name='prep_lumpy_evidence',
        ctx={'mem': 8, 'ncpus': 1},
        func="single_cell.utils.csvutils.rewrite_csv_file",
        args=(
            mgd.TempInputFile('lumpy_evidence.csv.gz'),
            mgd.OutputFile(lumpy_evidence, extensions=['.yaml']),
        ),
        kwargs={
            'dtypes': dtypes()['evidence']
        }
    )

    return workflow


def create_lumpy_workflow(
        config, normal_bam, tumour_cell_bams,
        lumpy_breakpoints_csv, lumpy_breakpoints_evidence,
        lumpy_breakpoints_bed
):
    ctx = {'docker_image': config['docker']['single_cell_pipeline']}
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=list(tumour_cell_bams.keys()),
    )

    workflow.subworkflow(
        name='normal_preprocess_lumpy',
        func='single_cell.workflows.lumpy.lumpy_preprocess_workflow',
        ctx={'docker_image': config['docker']['single_cell_pipeline']},
        args=(
            normal_bam,
            config,
            mgd.TempOutputFile('normal.discordants.sorted.bam'),
            mgd.TempOutputFile('normal.splitters.sorted.bam'),
            mgd.TempOutputFile('hist_normal_formatted.csv'),
            mgd.TempOutputFile('normal_mean_stdev.yaml')
        ),
    )

    workflow.subworkflow(
        name='tumour_preprocess_lumpy',
        func='single_cell.workflows.lumpy.lumpy_preprocess_workflow',
        ctx={'docker_image': config['docker']['single_cell_pipeline']},
        args=(
            mgd.InputFile('tumour_cells.bam','cell_id', extensions=['.bai'], fnames=tumour_cell_bams),
            config,
            mgd.TempOutputFile('tumour.discordants.sorted.bam'),
            mgd.TempOutputFile('tumour.splitters.sorted.bam'),
            mgd.TempOutputFile('hist_tumour_formatted.csv'),
            mgd.TempOutputFile('tumour_mean_stdev.yaml')
        ),
    )

    workflow.subworkflow(
        name='lumpy',
        ctx={'docker_image': config['docker']['single_cell_pipeline']},
        func="single_cell.workflows.lumpy.lumpy_calling_workflow",
        args=(
            config,
            mgd.TempInputFile('normal.discordants.sorted.bam'),
            mgd.TempInputFile('normal.splitters.sorted.bam'),
            mgd.TempInputFile('hist_normal_formatted.csv'),
            mgd.TempInputFile('normal_mean_stdev.yaml'),
            mgd.TempInputFile('tumour.discordants.sorted.bam'),
            mgd.TempInputFile('tumour.splitters.sorted.bam'),
            mgd.TempInputFile('hist_tumour_formatted.csv'),
            mgd.TempInputFile('tumour_mean_stdev.yaml'),
            mgd.OutputFile(lumpy_breakpoints_bed),
            mgd.OutputFile(lumpy_breakpoints_csv, extensions=['.yaml']),
            mgd.OutputFile(lumpy_breakpoints_evidence, extensions=['.yaml']),
        ),
    )

    return workflow
