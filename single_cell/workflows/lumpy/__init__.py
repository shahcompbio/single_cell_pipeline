import pypeliner
import pypeliner.managed as mgd


def lumpy_preprocess_cells(
        config, bam_files, merged_discordants, merged_splitters, hist_csv, mean_stdev_obj
):

    ctx = {'mem_retry_increment': 2, 'ncpus': 1,
           'docker_image': config['docker']['single_cell_pipeline']
           }

    lumpydocker = {'docker_image': config['docker']['lumpy']}

    histogram_settings = dict(
        N=10000, skip=100000, min_elements=1000, mads=10, X=4, read_length=101
    )
    histogram_settings.update(lumpydocker)

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=bam_files.keys(),
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


def create_lumpy_workflow(config, bam_files, normal_bam, lumpy_bed, lumpy_h5, tumour_id, normal_id):
    ctx = {'mem_retry_increment': 2, 'ncpus': 1,
           'docker_image': config['docker']['single_cell_pipeline']
           }

    lumpydocker = {'docker_image': config['docker']['lumpy']}

    histogram_settings = dict(
        N=10000, skip=100000, min_elements=1000, mads=10, X=4, read_length=101
    )
    histogram_settings.update(lumpydocker)

    print normal_bam

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('tumour_cell_id'),
        value=bam_files.keys(),
    )

    if isinstance(normal_bam, dict):
        workflow.setobj(
            obj=mgd.OutputChunks('normal_cell_id'),
            value=normal_bam.keys(),
        )
        workflow.set_filenames('normal_cells.bam', 'normal_cell_id', fnames=normal_bam)
        workflow.subworkflow(
            name='process_normal_cells',
            func='single_cell.workflows.lumpy.lumpy_preprocess_cells',
            args=(
                config,
                mgd.InputFile('normal_bam', 'normal_cell_id', fnames=normal_bam, extensions=['.bai']),
                mgd.TempOutputFile('normal.discordants.sorted.bam'),
                mgd.TempOutputFile('normal.splitters.sorted.bam'),
                mgd.TempOutputFile('hist_normal_formatted.csv'),
                mgd.TempOutputFile('normal_mean_stdev.yaml')
            ),
        )
    else:
        workflow.transform(
            name='process_normal',
            ctx={'mem': 8, 'ncpus': 1},
            func='single_cell.workflows.lumpy.tasks.process_bam',
            args=(
                mgd.InputFile(normal_bam, extensions=['.bai']),
                mgd.TempOutputFile('normal.discordants.sorted.bam'),
                mgd.TempOutputFile('normal.splitters.sorted.bam'),
                mgd.TempOutputFile('hist_normal.csv'),
                mgd.TempSpace("lumpy_normal_processing"),
            ),
            kwargs=histogram_settings,
        )
        workflow.transform(
            name='format_histo_normal',
            ctx={'mem': 8, 'ncpus': 1},
            func='single_cell.workflows.lumpy.merge_histograms.merge_histograms',
            args=(
                mgd.TempInputFile('hist_normal.csv'),
                mgd.TempOutputFile('hist_normal_formatted.csv'),
                mgd.TempOutputFile("normal_mean_stdev.yaml")
            ),
        )

    workflow.subworkflow(
        name='process_tumour_cells',
        func='single_cell.workflows.lumpy.lumpy_preprocess_cells',
        args=(
            config,
            mgd.InputFile('tumour_bam', 'tumour_cell_id', fnames=bam_files, extensions=['.bai']),
            mgd.TempOutputFile("tumour.discordants.sorted.bam"),
            mgd.TempOutputFile("tumour.splitters.sorted.bam"),
            mgd.TempOutputFile("tumour_hist.csv"),
            mgd.TempOutputFile('tumour_mean_stdev.yaml')
        ),
    )

    workflow.transform(
        name='run_lumpy',
        ctx={'mem': 8, 'ncpus': 1},
        func='single_cell.workflows.lumpy.tasks.run_lumpy',
        args=(
            mgd.TempInputFile("tumour.discordants.sorted.bam"),
            mgd.TempInputFile("tumour.splitters.sorted.bam"),
            mgd.TempInputFile("tumour_hist.csv"),
            mgd.TempInputFile('tumour_mean_stdev.yaml'),
            tumour_id,
            mgd.TempInputFile('normal.discordants.sorted.bam'),
            mgd.TempInputFile('normal.splitters.sorted.bam'),
            mgd.TempInputFile('hist_normal_formatted.csv'),
            mgd.TempInputFile('normal_mean_stdev.yaml'),
            normal_id,
            mgd.OutputFile(lumpy_bed),
            mgd.TempSpace("lumpy_temp"),
        ),
        kwargs=lumpydocker,
    )

    workflow.transform(
        name='parse_lumpy',
        ctx={'mem': 8, 'ncpus': 1},
        func='single_cell.workflows.lumpy.parse_lumpy_to_h5.parse_lumpy',
        args=(
            mgd.InputFile(lumpy_bed),
            mgd.OutputFile(lumpy_h5),
        ),
    )

    return workflow
