'''
Created on Jul 6, 2017

@author: dgrewal
'''

import pypeliner
import pypeliner.managed as mgd


def bam_metrics_workflow(
        bam_filename,
        summary_fastq_screen_count_per_cell,
        alignment_metrics,
        gc_metrics,
        markdups_metrics_percell,
        flagstat_metrics_percell,
        wgs_metrics_percell,
        gc_metrics_percell,
        gc_metrics_summary_percell,
        gc_metrics_pdf_percell,
        insert_metrics_percell,
        insert_metrics_pdf_percell,
        ref_genome,
        sample_info,
        config,
        cell_ids
):
    markdups_metrics_percell = dict([(cellid, markdups_metrics_percell[cellid])
                                     for cellid in cell_ids])

    flagstat_metrics_percell = dict([(cellid, flagstat_metrics_percell[cellid])
                                     for cellid in cell_ids])

    wgs_metrics_percell = dict([(cellid, wgs_metrics_percell[cellid])
                                for cellid in cell_ids])

    gc_metrics_percell = dict([(cellid, gc_metrics_percell[cellid])
                               for cellid in cell_ids])

    gc_metrics_summary_percell = dict([(cellid, gc_metrics_summary_percell[cellid])
                                       for cellid in cell_ids])

    gc_metrics_pdf_percell = dict([(cellid, gc_metrics_pdf_percell[cellid])
                                   for cellid in cell_ids])

    insert_metrics_percell = dict([(cellid, insert_metrics_percell[cellid])
                                   for cellid in cell_ids])

    insert_metrics_pdf_percell = dict([(cellid, insert_metrics_pdf_percell[cellid])
                                       for cellid in cell_ids])

    baseimage = config['docker']['single_cell_pipeline']
    workflow = pypeliner.workflow.Workflow(ctx={'docker_image': baseimage})

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cell_ids,
    )

    workflow.transform(
        name='get_duplication_metrics',
        axes=('cell_id',),
        func="single_cell.utils.picardutils.bam_markdups",
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            mgd.TempOutputFile("temp_markdup_bam.bam", 'cell_id'),
            mgd.OutputFile('markdups_metrics', 'cell_id', fnames=markdups_metrics_percell),
            mgd.TempSpace('tempdir_markdups', 'cell_id'),
        ),
        kwargs={'docker_image': config['docker']['picard']}
    )

    workflow.transform(
        name='get_flagstat_metrics',
        axes=('cell_id',),
        func="single_cell.utils.bamutils.bam_flagstat",
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            mgd.OutputFile('flagstat_metrics_percell', 'cell_id', fnames=flagstat_metrics_percell),
        ),
        kwargs={'docker_image': config['docker']['samtools']}
    )

    workflow.transform(
        name='bam_collect_wgs_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        func="single_cell.utils.picardutils.bam_collect_wgs_metrics",
        axes=('cell_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            ref_genome,
            mgd.OutputFile('wgs_metrics_percell', 'cell_id', fnames=wgs_metrics_percell),
            config['picard_wgs_params'],
            mgd.TempSpace('wgs_tempdir', 'cell_id'),
        ),
        kwargs={
            'docker_image': config['docker']['picard']
        }
    )

    workflow.transform(
        name='bam_collect_gc_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        func="single_cell.utils.picardutils.bam_collect_gc_metrics",
        axes=('cell_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            ref_genome,
            mgd.OutputFile('gc_metrics_percell', 'cell_id', fnames=gc_metrics_percell),
            mgd.OutputFile('gc_metrics_summary_percell', 'cell_id', fnames=gc_metrics_summary_percell),
            mgd.OutputFile('gc_metrics_pdf_percell', 'cell_id', fnames=gc_metrics_pdf_percell),
            mgd.TempSpace('gc_tempdir', 'cell_id'),
        ),
        kwargs={
            'docker_image': config['docker']['picard']
        }
    )

    workflow.transform(
        name='bam_collect_insert_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        func="single_cell.utils.picardutils.bam_collect_insert_metrics",
        axes=('cell_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            mgd.InputFile('flagstat_metrics_percell', 'cell_id', fnames=flagstat_metrics_percell),
            mgd.OutputFile('insert_metrics_percell', 'cell_id', fnames=insert_metrics_percell),
            mgd.OutputFile('insert_metrics_pdf_percell', 'cell_id', fnames=insert_metrics_pdf_percell),
            mgd.TempSpace('insert_tempdir', 'cell_id'),
        ),
        kwargs={
            'docker_image': config['docker']['picard']
        }
    )

    workflow.transform(
        name="collect_gc_metrics",
        func="single_cell.workflows.align.tasks.collect_gc",
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        args=(
            mgd.InputFile('gc_metrics_percell', 'cell_id', axes_origin=[], fnames=gc_metrics_percell),
            mgd.OutputFile(gc_metrics, extensions=['.yaml']),
            mgd.TempSpace("temp_gc")
        ),
    )

    workflow.transform(
        name='collect_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        func="single_cell.workflows.align.tasks.collect_metrics",
        args=(
            mgd.InputFile('flagstat_metrics', 'cell_id', axes_origin=[], fnames=flagstat_metrics_percell),
            mgd.InputFile('markdups_metrics', 'cell_id', axes_origin=[], fnames=markdups_metrics_percell),
            mgd.InputFile('insert_metrics_percell', 'cell_id', axes_origin=[], fnames=insert_metrics_percell),
            mgd.InputFile('wgs_metrics_percell', 'cell_id', axes_origin=[], fnames=wgs_metrics_percell),
            mgd.TempSpace("tempdir_collect_metrics"),
            mgd.TempOutputFile("alignment_metrics.csv.gz", extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='annotate_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        func="single_cell.utils.csvutils.annotate_csv",
        args=(
            mgd.TempInputFile("alignment_metrics.csv.gz", extensions=['.yaml']),
            sample_info,
            mgd.TempOutputFile('alignment_metrics_annotated.csv.gz', extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='add_fastqscreen_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        func="single_cell.utils.csvutils.merge_csv",
        args=(
            [
                mgd.TempInputFile("alignment_metrics_annotated.csv.gz", extensions=['.yaml']),
                mgd.InputFile(summary_fastq_screen_count_per_cell),
            ],
            mgd.OutputFile(alignment_metrics, extensions=['.yaml']),
            'outer',
            ['cell_id'],
        ),
    )

    return workflow


def create_alignment_workflow(
        fastq_1_filename,
        fastq_2_filename,
        bam_filename,
        alignment_metrics,
        gc_metrics,
        detailed_fastqscreen_metrics,
        plot_metrics,
        ref_genome,
        config,
        triminfo,
        centerinfo,
        sample_info,
        cell_ids,
        metrics_tar,
        library_id,
        realign=False
):
    baseimage = config['docker']['single_cell_pipeline']

    bam_filename = dict([(cellid, bam_filename[cellid])
                         for cellid in cell_ids])

    chromosomes = config["chromosomes"]

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('chrom'),
        value=chromosomes,
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('sampleinfo', 'cell_id', axes_origin=[]),
        value=sample_info)

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id', 'lane'),
        value=fastq_1_filename.keys(),
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('trim', 'cell_id', 'lane', axes_origin=[]),
        value=triminfo)

    workflow.setobj(
        obj=mgd.TempOutputObj('center', 'cell_id', 'lane', axes_origin=[]),
        value=centerinfo)

    workflow.transform(
        name='run_fastq_screen',
        ctx={'mem': 7, 'ncpus': 1, 'docker_image': baseimage},
        axes=('cell_id', 'lane',),
        func="single_cell.workflows.align.fastqscreen.organism_filter",
        args=(
            mgd.InputFile('fastq_1', 'cell_id', 'lane', fnames=fastq_1_filename),
            mgd.InputFile('fastq_2', 'cell_id', 'lane', fnames=fastq_2_filename),
            mgd.TempOutputFile('fastq_r1_matching_reads.fastq.gz', 'cell_id', 'lane'),
            mgd.TempOutputFile('fastq_r2_matching_reads.fastq.gz', 'cell_id', 'lane'),
            mgd.TempOutputFile('organism_detailed_count_per_lane.csv', 'cell_id', 'lane'),
            mgd.TempOutputFile('organism_summary_count_per_lane.csv', 'cell_id', 'lane'),
            mgd.TempSpace("tempdir_organism_filter", 'cell_id', 'lane'),
            mgd.InputInstance('cell_id'),
            config['fastq_screen_params'],
            config['ref_type']
        ),
        kwargs={
            'docker_image': config['docker']['fastq_screen'],
            'filter_contaminated_reads': config['fastq_screen_params']['filter_contaminated_reads']
        }
    )

    workflow.transform(
        name='merge_fastq_screen_metrics',
        ctx={'mem': 7, 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.align.fastqscreen.merge_fastq_screen_counts",
        args=(
            mgd.TempInputFile('organism_detailed_count_per_lane.csv', 'cell_id', 'lane'),
            mgd.TempInputFile('organism_summary_count_per_lane.csv', 'cell_id', 'lane'),
            mgd.OutputFile(detailed_fastqscreen_metrics, extensions=['.yaml']),
            mgd.TempOutputFile('organism_summary_count_per_cell.csv'),
        )
    )

    workflow.transform(
        name='align_reads',
        ctx={'mem': 7, 'ncpus': 1, 'docker_image': baseimage},
        axes=('cell_id', 'lane',),
        func="single_cell.workflows.align.tasks.align_pe",
        args=(
            mgd.TempInputFile('fastq_r1_matching_reads.fastq.gz', 'cell_id', 'lane'),
            mgd.TempInputFile('fastq_r2_matching_reads.fastq.gz', 'cell_id', 'lane'),
            mgd.TempOutputFile(
                'aligned_per_cell_per_lane.sorted.bam', 'cell_id', 'lane'),
            mgd.TempOutputFile('fastqc_reports.tar.gz', 'cell_id', 'lane'),
            mgd.TempOutputFile('flagstat_metrics.txt', 'cell_id', 'lane'),
            mgd.TempSpace('alignment_temp', 'cell_id', 'lane'),
            ref_genome,
            mgd.TempInputObj('trim', 'cell_id', 'lane'),
            mgd.TempInputObj('center', 'cell_id', 'lane'),
            mgd.TempInputObj('sampleinfo', 'cell_id'),
            mgd.InputInstance('cell_id'),
            mgd.InputInstance('lane'),
            library_id,
            config['aligner'],
            config['docker'],
            config['adapter'],
            config['adapter2'],
            config['fastq_screen_params'],
        )
    )

    workflow.transform(
        name='merge_bams',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.align.tasks.merge_bams",
        axes=('cell_id',),
        args=(
            mgd.TempInputFile(
                'aligned_per_cell_per_lane.sorted.bam',
                'cell_id',
                'lane'),
            mgd.TempOutputFile('merged_lanes.bam', 'cell_id'),
            mgd.TempOutputFile('merged_lanes.bam.bai', 'cell_id'),
            config['docker']
        )
    )

    if realign:
        workflow.transform(
            name='realignment',
            axes=('chrom',),
            ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
            func="single_cell.workflows.align.tasks.realign",
            args=(
                mgd.TempInputFile('merged_lanes.bam', 'cell_id'),
                mgd.TempInputFile('merged_lanes.bam.bai', 'cell_id'),
                mgd.TempOutputFile('realigned.bam', 'chrom', 'cell_id'),
                mgd.TempSpace('realignment_temp', 'chrom'),
                config,
                mgd.InputInstance('chrom')
            )
        )

        workflow.transform(
            name='merge_realignment',
            ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
            axes=('cell_id',),
            func="single_cell.workflows.align.tasks.merge_realignment",
            args=(
                mgd.TempInputFile('realigned.bam', 'chrom', 'cell_id'),
                mgd.TempOutputFile('merged_realign.bam', 'cell_id'),
                config,
                mgd.InputInstance('cell_id')
            )
        )

    final_bam = mgd.TempInputFile('merged_lanes.bam', 'cell_id')
    if realign:
        final_bam = mgd.TempInputFile('merged_realign.bam', 'cell_id')

    workflow.transform(
        name='postprocess_bam',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        axes=('cell_id',),
        func="single_cell.workflows.align.tasks.postprocess_bam",
        args=(
            final_bam,
            mgd.OutputFile('sorted_markdups', 'cell_id', fnames=bam_filename, extensions=['.bai']),
            mgd.TempSpace('tempdir', 'cell_id'),
            config['docker'],
        ),
    )

    workflow.subworkflow(
        name='metrics_subworkflow',
        func="single_cell.workflows.align.bam_metrics_workflow",
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename, extensions=['.bai']),
            mgd.TempInputFile('organism_summary_count_per_cell.csv'),
            mgd.OutputFile(alignment_metrics, extensions=['.yaml']),
            mgd.OutputFile(gc_metrics, extensions=['.yaml']),
            mgd.TempOutputFile('markdups_metrics.txt', 'cell_id', axes_origin=[]),
            mgd.TempOutputFile('flagstat_metrics.txt', 'cell_id', axes_origin=[]),
            mgd.TempOutputFile('wgs_metrics.txt', 'cell_id', axes_origin=[]),
            mgd.TempOutputFile('gc_metrics.txt', 'cell_id', axes_origin=[]),
            mgd.TempOutputFile('gc_metrics_summary.txt', 'cell_id', axes_origin=[]),
            mgd.TempOutputFile('gc_metrics.pdf', 'cell_id', axes_origin=[]),
            mgd.TempOutputFile('insert_metrics.txt', 'cell_id', axes_origin=[]),
            mgd.TempOutputFile('insert_metrics.pdf', 'cell_id', axes_origin=[]),
            ref_genome,
            sample_info,
            config,
            cell_ids
        )
    )

    workflow.transform(
        name='plot_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.align.tasks.plot_metrics",
        args=(
            mgd.InputFile(alignment_metrics, extensions=['.yaml']),
            mgd.OutputFile(plot_metrics),
            'QC pipeline metrics',
            mgd.InputFile(gc_metrics, extensions=['.yaml']),
            config['gc_windows'],
        )
    )

    workflow.transform(
        name='tar_all_files',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.utils.helpers.tar_files",
        args=(
            [
                mgd.TempInputFile('fastqc_reports.tar.gz', 'cell_id', 'lane'),
                mgd.TempInputFile('flagstat_metrics.txt', 'cell_id', 'lane'),
                mgd.TempInputFile('markdups_metrics.txt', 'cell_id'),
                mgd.TempInputFile('flagstat_metrics.txt', 'cell_id'),
                mgd.TempInputFile('wgs_metrics.txt', 'cell_id'),
                mgd.TempInputFile('gc_metrics.txt', 'cell_id'),
                mgd.TempInputFile('gc_metrics_summary.txt', 'cell_id'),
                mgd.TempInputFile('gc_metrics.pdf', 'cell_id'),
                mgd.TempInputFile('insert_metrics.txt', 'cell_id'),
                mgd.TempInputFile('insert_metrics.pdf', 'cell_id'),
            ],
            mgd.OutputFile(metrics_tar),
            mgd.TempSpace("merge_metrics_tar")
        )
    )

    return workflow
