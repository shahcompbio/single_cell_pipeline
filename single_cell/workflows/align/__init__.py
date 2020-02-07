'''
Created on Jul 6, 2017

@author: dgrewal
'''

import pypeliner
import pypeliner.managed as mgd
from single_cell.workflows.align.dtypes import dtypes

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
        name='get_duplication_wgs_flagstat_metrics',
        axes=('cell_id',),
        func="single_cell.workflows.align.tasks.picard_wgs_dup",
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            mgd.TempOutputFile("temp_markdup_bam.bam", 'cell_id'),
            mgd.OutputFile('markdups_metrics', 'cell_id', fnames=markdups_metrics_percell),
            mgd.TempSpace('tempdir_markdups', 'cell_id'),
            ref_genome,
            mgd.OutputFile('wgs_metrics_percell', 'cell_id', fnames=wgs_metrics_percell),
            config['picard_wgs_params'],
        ),
        kwargs={
            'picard_docker': config['docker']['picard'],
        }
    )

    workflow.transform(
        name='bam_collect_gc_insert_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        func="single_cell.workflows.align.tasks.picard_insert_gc_flagstat",
        axes=('cell_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            ref_genome,
            mgd.OutputFile('gc_metrics_percell', 'cell_id', fnames=gc_metrics_percell),
            mgd.OutputFile('gc_metrics_summary_percell', 'cell_id', fnames=gc_metrics_summary_percell),
            mgd.OutputFile('gc_metrics_pdf_percell', 'cell_id', fnames=gc_metrics_pdf_percell),
            mgd.TempSpace('gc_tempdir', 'cell_id'),
            mgd.OutputFile('flagstat_metrics_percell', 'cell_id', fnames=flagstat_metrics_percell),
            mgd.OutputFile('insert_metrics_percell', 'cell_id', fnames=insert_metrics_percell),
            mgd.OutputFile('insert_metrics_pdf_percell', 'cell_id', fnames=insert_metrics_pdf_percell),
        ),
        kwargs={
            'picard_docker': config['docker']['picard'],
            'samtools_docker': config['docker']['samtools']
        }
    )

    workflow.transform(
        name="collect_gc_metrics",
        func="single_cell.workflows.align.tasks.collect_gc",
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        args=(
            mgd.InputFile('gc_metrics_percell.csv.gz', 'cell_id', axes_origin=[], fnames=gc_metrics_percell),
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
        kwargs={'annotation_dtypes': dtypes()['metrics']}
    )

    workflow.transform(
        name='add_fastqscreen_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        func="single_cell.utils.csvutils.merge_csv",
        args=(
            [
                mgd.TempInputFile("alignment_metrics_annotated.csv.gz", extensions=['.yaml']),
                mgd.InputFile(summary_fastq_screen_count_per_cell, extensions=['.yaml']),
            ],
            mgd.OutputFile(alignment_metrics, extensions=['.yaml']),
            'outer',
            ['cell_id'],
        ),
        # kwargs={'dtypes': dtypes()['metrics']}
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
        laneinfo,
        sample_info,
        cell_ids,
        metrics_tar,
        library_id,
):

    baseimage = config['docker']['single_cell_pipeline']

    ctx = {'mem': 7, 'ncpus': 1, 'docker_image': baseimage, 'mem_retry_factor': 1}

    bam_filename = dict([(cellid, bam_filename[cellid])
                         for cellid in cell_ids])

    chromosomes = config["chromosomes"]

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('chrom'),
        value=chromosomes,
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('sampleinfo', 'cell_id', axes_origin=[]),
        value=sample_info)

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id', 'lane'),
        value=list(fastq_1_filename.keys()),
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('laneinfo', 'cell_id', 'lane', axes_origin=[]),
        value=laneinfo)

    workflow.transform(
        name='align_reads',
        axes=('cell_id',),
        func="single_cell.workflows.align.align_tasks.align_lanes",
        args=(
            mgd.InputFile('fastq_1', 'cell_id', 'lane', fnames=fastq_1_filename, axes_origin=[]),
            mgd.InputFile('fastq_2', 'cell_id', 'lane', fnames=fastq_2_filename, axes_origin=[]),
            mgd.OutputFile('sorted_markdups', 'cell_id', fnames=bam_filename, extensions=['.bai']),
            mgd.TempOutputFile('fastqc_reports.tar.gz', 'cell_id'),
            mgd.TempSpace('alignment_temp', 'cell_id'),
            ref_genome,
            mgd.TempInputObj('laneinfo', 'cell_id', 'lane', axes_origin=[]),
            mgd.TempInputObj('sampleinfo', 'cell_id'),
            mgd.InputInstance('cell_id'),
            library_id,
            config['aligner'],
            config['docker'],
            config['adapter'],
            config['adapter2'],
            mgd.TempOutputFile('organism_detailed_count_per_cell.csv.gz', 'cell_id'),
            mgd.TempOutputFile('organism_summary_count_per_cell.csv.gz', 'cell_id'),
            config['fastq_screen_params'],
        )
    )

    workflow.transform(
        name='merge_fastq_screen_metrics',
        func="single_cell.workflows.align.fastqscreen.merge_fastq_screen_counts",
        args=(
            mgd.TempInputFile('organism_detailed_count_per_cell.csv.gz', 'cell_id'),
            mgd.TempInputFile('organism_summary_count_per_cell.csv.gz', 'cell_id'),
            mgd.OutputFile(detailed_fastqscreen_metrics, extensions=['.yaml']),
            mgd.TempOutputFile('organism_summary_count_per_cell.csv.gz', extensions=['.yaml']),
        )
    )

    workflow.subworkflow(
        name='metrics_subworkflow',
        func="single_cell.workflows.align.bam_metrics_workflow",
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename, extensions=['.bai']),
            mgd.TempInputFile('organism_summary_count_per_cell.csv.gz', extensions=['.yaml']),
            mgd.TempOutputFile('alignment_metrics.csv.gz', extensions=['.yaml']),
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
        name='add_contamination_status',
        ctx={'mem': config['memory']['med']},
        func="single_cell.workflows.align.tasks.add_contamination_status",
        args=(
            mgd.TempInputFile('alignment_metrics.csv.gz', extensions=['.yaml']),
            mgd.OutputFile(alignment_metrics, extensions=['.yaml']),
        ),
        kwargs={
            'reference': config['ref_type'],
            'strict_validation': config['fastq_screen_params']['strict_validation']
        }
    )

    workflow.transform(
        name='plot_metrics',
        ctx={'mem': config['memory']['med']},
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
        ctx={'mem': config['memory']['med']},
        func="single_cell.utils.helpers.tar_files",
        args=(
            [
                mgd.TempInputFile('fastqc_reports.tar.gz', 'cell_id'),
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
