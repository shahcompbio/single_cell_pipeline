'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
import single_cell
from single_cell.utils import helpers

def create_alignment_metrics_workflow(
        bam_filename,
        alignment_metrics,
        gc_metrics,
        plot_metrics,
        ref_genome,
        config,
        args,
        sample_info,
        cell_ids):

    baseimage = config['docker']['single_cell_pipeline']
    out_dir = args['out_dir']

    merge_metrics = os.path.join(out_dir, 'metrics')

    bam_filename = dict([(cellid, bam_filename[cellid])
                         for cellid in cell_ids])

    chromosomes = config["chromosomes"]

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('chrom'),
        value=chromosomes,
    )

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cell_ids,
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('sampleinfo', 'cell_id', axes_origin=[]),
        value=sample_info)


    markdups_metrics = os.path.join(
        merge_metrics,
        'markdups_metrics',
        '{cell_id}.markdups_metrics.txt')
    flagstat_metrics = os.path.join(
        merge_metrics,
        'flagstat_metrics',
        '{cell_id}.flagstat_metrics.txt')
    workflow.transform(
        name='postprocess_bam_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        axes=('cell_id',),
        func="single_cell.workflows.alignment_metrics.tasks.get_postprocess_metrics",
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename, extensions=['.bai']),
            mgd.TempSpace('tempdir', 'cell_id'),
            config['docker'],
            mgd.OutputFile(markdups_metrics, 'cell_id'),
            mgd.OutputFile(flagstat_metrics, 'cell_id'),
        ),
    )


    wgs_metrics_filename = os.path.join(
        merge_metrics, 'wgs_metrics', '{cell_id}.wgs_metrics.txt')
    workflow.transform(
        name='bam_collect_wgs_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.alignment_metrics.tasks.bam_collect_wgs_metrics",
        axes=('cell_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            ref_genome,
            mgd.OutputFile(wgs_metrics_filename, 'cell_id'),
            config['docker'],
            config['picard_wgs_params'],
            mgd.TempSpace('wgs_tempdir', 'cell_id'),
        ),
    )

    gc_metrics_filename = os.path.join(
        merge_metrics,
        'gc_metrics',
        '{cell_id}.gc_metrics.txt')
    gc_summary_filename = os.path.join(
        merge_metrics,
        'gc_metrics',
        '{cell_id}.gc_metrics.summ.txt')
    gc_chart_filename = os.path.join(
        merge_metrics,
        'gc_metrics',
        '{cell_id}.gc_metrics.pdf')
    workflow.transform(
        name='bam_collect_gc_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.alignment_metrics.tasks.bam_collect_gc_metrics",
        axes=('cell_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            ref_genome,
            mgd.OutputFile(gc_metrics_filename, 'cell_id'),
            mgd.OutputFile(gc_summary_filename, 'cell_id'),
            mgd.OutputFile(gc_chart_filename, 'cell_id'),
            mgd.TempSpace('gc_tempdir', 'cell_id'),
            config['docker'],
        ),
    )

    insert_metrics_filename = os.path.join(
        merge_metrics,
        'insert_metrics',
        '{cell_id}.insert_metrics.txt')
    insert_histogram_filename = os.path.join(
        merge_metrics,
        'insert_metrics',
        '{cell_id}.insert_metrics.pdf')
    workflow.transform(
        name='bam_collect_insert_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.alignment_metrics.tasks.bam_collect_insert_metrics",
        axes=('cell_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            mgd.InputFile(flagstat_metrics, 'cell_id'),
            mgd.OutputFile(insert_metrics_filename, 'cell_id'),
            mgd.OutputFile(insert_histogram_filename, 'cell_id'),
            mgd.TempSpace('insert_tempdir', 'cell_id'),
            config['docker'],
        ),
    )

    workflow.transform(
        name="collect_gc_metrics",
        func="single_cell.workflows.alignment_metrics.tasks.collect_gc",
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        args=(
            mgd.InputFile(gc_metrics_filename, 'cell_id', axes_origin=[]),
            mgd.TempOutputFile('gc_metrics.csv.gz', extensions=['.yaml']),
            mgd.TempSpace("temp_gc")
        ),
    )

    workflow.transform(
        name='collect_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.alignment_metrics.tasks.collect_metrics",
        args=(
            mgd.InputFile(flagstat_metrics, 'cell_id', axes_origin=[]),
            mgd.InputFile(markdups_metrics, 'cell_id', axes_origin=[]),
            mgd.InputFile(insert_metrics_filename, 'cell_id', axes_origin=[]),
            mgd.InputFile(wgs_metrics_filename, 'cell_id', axes_origin=[]),
            mgd.TempSpace("tempdir_collect_metrics"),
            mgd.TempOutputFile("alignment_metrics.csv.gz", extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='annotate_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.alignment_metrics.tasks.annotate_metrics",
        args=(
            mgd.TempInputFile("alignment_metrics.csv.gz", extensions=['.yaml']),
            sample_info,
            mgd.TempOutputFile('alignment_metrics_annotated.csv.gz', extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='plot_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.alignment_metrics.tasks.plot_metrics",
        args=(
            mgd.TempInputFile('alignment_metrics_annotated.csv.gz', extensions=['.yaml']),
            mgd.OutputFile(plot_metrics),
            'QC pipeline metrics',
            mgd.TempInputFile('gc_metrics.csv.gz', extensions=['.yaml']),
            config['gc_windows'],
        )
    )

    workflow.transform(
        name='finalize_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        axes=('multiplier',),
        func="single_cell.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile('alignment_metrics_annotated.csv.gz', extensions=['.yaml']),
            mgd.OutputFile(alignment_metrics, extensions=['.yaml']),
        ),
    )


    workflow.transform(
        name='finalize_gc_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        axes=('multiplier',),
        func="single_cell.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile('gc_metrics.csv.gz', extensions=['.yaml']),
            mgd.OutputFile(gc_metrics, extensions=['.yaml']),
        ),
    )


    return workflow
