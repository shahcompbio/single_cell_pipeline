import pypeliner
import pypeliner.managed as mgd


def create_lumpy_workflow(config, bam_files, normal_bam, lumpy_bed, lumpy_h5, tumour_id, normal_id):
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
        obj=mgd.OutputChunks('sample_id'),
        value=bam_files.keys(),
    )

    workflow.transform(
        name='process_normal',
        ctx={'mem': 8, 'ncpus': 1},
        func='single_cell.workflows.lumpy.tasks.process_bam',
        args=(
            mgd.InputFile(normal_bam),
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
        ret=mgd.TempOutputObj('normal_mean_stdev'),
        func='single_cell.workflows.lumpy.merge_histograms.merge_histograms',
        args=(
            mgd.TempInputFile('hist_normal.csv'),
            mgd.TempOutputFile('hist_normal_formatted.csv'),
        ),
    )

    workflow.transform(
        name='process_tumour_cells',
        axes=('sample_id',),
        ctx={'mem': 8, 'ncpus': 1},
        func='single_cell.workflows.lumpy.tasks.process_bam',
        args=(
            mgd.InputFile('tumour_bam', 'sample_id', fnames=bam_files),
            mgd.TempOutputFile('tumour.discordants.sorted.bam', 'sample_id'),
            mgd.TempOutputFile('tumour.splitters.sorted.bam', 'sample_id'),
            mgd.TempOutputFile('hist.csv', 'sample_id'),
            mgd.TempSpace("lumpy_tumour_processing", "sample_id"),
        ),
        kwargs=dict(tag=mgd.InputInstance('sample_id'), **histogram_settings),
    )

    workflow.transform(
        name='merge_disc',
        ctx={'mem': 8, 'ncpus': 1},
        func='single_cell.workflows.lumpy.tasks.merge_bams',
        args=(
            mgd.TempInputFile('tumour.discordants.sorted.bam', 'sample_id'),
            mgd.TempOutputFile("merged_discordants.sorted.bam")
        ),
        kwargs=lumpydocker,
    )

    workflow.transform(
        name='merge_split',
        ctx={'mem': 8, 'ncpus': 1},
        func='single_cell.workflows.lumpy.tasks.merge_bams',
        args=(
            mgd.TempInputFile('tumour.splitters.sorted.bam', 'sample_id'),
            mgd.TempOutputFile("merged_splitters.sorted.bam")
        ),
        kwargs=lumpydocker,
    )

    workflow.transform(
        name='merge_histo',
        ctx={'mem': 8, 'ncpus': 1},
        ret=mgd.TempOutputObj('tumour_mean_stdev'),
        func='single_cell.workflows.lumpy.merge_histograms.merge_histograms',
        args=(
            mgd.TempInputFile('hist.csv', 'sample_id'),
            mgd.TempOutputFile("hist.csv")
        ),
    )


    workflow.transform(
        name='run_lumpy',
        ctx={'mem': 8, 'ncpus': 1},
        func='single_cell.workflows.lumpy.tasks.run_lumpy',
        args=(
            mgd.TempInputFile("merged_discordants.sorted.bam"),
            mgd.TempInputFile("merged_splitters.sorted.bam"),
            mgd.TempInputFile("hist.csv"),
            mgd.TempInputObj('tumour_mean_stdev'),
            tumour_id,
            mgd.TempInputFile('normal.discordants.sorted.bam'),
            mgd.TempInputFile('normal.splitters.sorted.bam'),
            mgd.TempInputFile('hist_normal_formatted.csv'),
            mgd.TempInputObj('normal_mean_stdev'),
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
