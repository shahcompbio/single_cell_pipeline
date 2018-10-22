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
        bai_filename,
        alignment_metrics,
        plot_metrics,
        ref_genome,
        config,
        args,
        sample_info,
        cell_ids):

    out_dir = args['out_dir']

    merge_metrics = os.path.join(out_dir, 'metrics')

    bam_filename = dict([(cellid, bam_filename[cellid])
                         for cellid in cell_ids])

    bai_filename = dict([(cellid, bai_filename[cellid])
                         for cellid in cell_ids])

    chromosomes = config["chromosomes"]

    ctx = {'mem_retry_increment': 2, 'ncpus': 1}
    docker_ctx = helpers.get_container_ctx(config['containers'], 'single_cell_pipeline')
    ctx.update(docker_ctx)

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
        name='postprocess_bam',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        axes=('cell_id',),
        func="single_cell.workflows.alignment_metrics.tasks.get_postprocess_metrics",
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bai_filename),
            mgd.TempSpace('tempdir', 'cell_id'),
            config,
            mgd.OutputFile(markdups_metrics, 'cell_id'),
            mgd.OutputFile(flagstat_metrics, 'cell_id'),
        ),
    )


    wgs_metrics_filename = os.path.join(
        merge_metrics, 'wgs_metrics', '{cell_id}.wgs_metrics.txt')
    workflow.transform(
        name='bam_collect_wgs_metrics',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        func="single_cell.workflows.alignment_metrics.tasks.bam_collect_wgs_metrics",
        axes=('cell_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            ref_genome,
            mgd.OutputFile(wgs_metrics_filename, 'cell_id'),
            config,
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
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        func="single_cell.workflows.alignment_metrics.tasks.bam_collect_gc_metrics",
        axes=('cell_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            ref_genome,
            mgd.OutputFile(gc_metrics_filename, 'cell_id'),
            mgd.OutputFile(gc_summary_filename, 'cell_id'),
            mgd.OutputFile(gc_chart_filename, 'cell_id'),
            mgd.TempSpace('gc_tempdir', 'cell_id'),
            config,
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
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        func="single_cell.workflows.alignment_metrics.tasks.bam_collect_insert_metrics",
        axes=('cell_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            mgd.InputFile(flagstat_metrics, 'cell_id'),
            mgd.OutputFile(insert_metrics_filename, 'cell_id'),
            mgd.OutputFile(insert_histogram_filename, 'cell_id'),
            mgd.TempSpace('insert_tempdir', 'cell_id'),
            config,
        ),
    )

    workflow.transform(
        name="collect_gc_metrics",
        func="single_cell.workflows.alignment_metrics.tasks.collect_gc",
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        args=(
            mgd.InputFile(gc_metrics_filename, 'cell_id', axes_origin=[]),
            mgd.TempOutputFile("gc_metrics.h5"),
            mgd.TempSpace("temp_gc")
        )
    )

    workflow.transform(
        name='collect_metrics',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        func="single_cell.workflows.alignment_metrics.tasks.collect_metrics",
        args=(
            mgd.InputFile(flagstat_metrics, 'cell_id', axes_origin=[]),
            mgd.InputFile(markdups_metrics, 'cell_id', axes_origin=[]),
            mgd.InputFile(insert_metrics_filename, 'cell_id', axes_origin=[]),
            mgd.InputFile(wgs_metrics_filename, 'cell_id', axes_origin=[]),
            mgd.TempSpace("tempdir_collect_metrics"),
            mgd.TempOutputFile("alignment_metrics.h5"),
        ),
    )

    workflow.transform(
        name='annotate_metrics',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        func="single_cell.workflows.alignment_metrics.tasks.annotate_metrics",
        args=(
            mgd.TempInputFile("alignment_metrics.h5"),
            sample_info,
            mgd.TempOutputFile("alignment_metrics_annotated.h5"),
        )
    )

    workflow.transform(
        name='plot_metrics',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        func="single_cell.workflows.alignment_metrics.tasks.plot_metrics",
        args=(
            mgd.TempInputFile("alignment_metrics_annotated.h5"),
            mgd.OutputFile(plot_metrics),
            'QC pipeline metrics',
            mgd.TempInputFile("gc_metrics.h5"),
            config['gc_windows'],
        )
    )

    workflow.transform(
        name='concatenate_all_hdf_tables',
        ctx=dict(mem=config['memory']['low'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        func="single_cell.utils.hdfutils.concat_hdf_tables",
        args=(
            [mgd.TempInputFile("alignment_metrics_annotated.h5"),
             mgd.TempInputFile("gc_metrics.h5"),
             ],
            mgd.TempOutputFile("alignment_precast.h5"),
        ),
    )

    workflow.transform(
        name='cast_h5',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        func="single_cell.utils.hdfutils.cast_h5_file",
        args=(
            mgd.TempInputFile("alignment_precast.h5"),
            mgd.OutputFile(alignment_metrics),
        )
    )


    return workflow
