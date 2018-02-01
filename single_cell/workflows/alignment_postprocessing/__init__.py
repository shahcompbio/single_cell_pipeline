'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
import tasks

def create_bam_post_workflow(
    bam,
    bam_filename,
    bam_index_filename,
    ref_genome,
    alignment_metrics,
    gc_metrics,
    sample_ids,
    config,
    out_dir):
 
    metrics_dir = os.path.join(out_dir, 'metrics')
    markdups_metrics_filename = os.path.join(metrics_dir, 'markdups_metrics', '{sample_id}.markdups_metrics.txt')
    flagstat_metrics_filename = os.path.join(metrics_dir, 'flagstat_metrics', '{sample_id}.flagstat_metrics.txt')
    wgs_metrics_filename = os.path.join(metrics_dir, 'wgs_metrics', '{sample_id}.wgs_metrics.txt')
    gc_metrics_filename = os.path.join(metrics_dir, 'gc_metrics', '{sample_id}.gc_metrics.txt')
    gc_summary_filename = os.path.join(metrics_dir, 'gc_metrics', '{sample_id}.gc_metrics.summ.txt')
    gc_chart_filename = os.path.join(metrics_dir, 'gc_metrics', '{sample_id}.gc_metrics.pdf')
    insert_metrics_filename = os.path.join(metrics_dir, 'insert_metrics', '{sample_id}.insert_metrics.txt')
    insert_histogram_filename = os.path.join(metrics_dir, 'insert_metrics', '{sample_id}.insert_metrics.pdf')



    bam_filename = dict([(sampid, bam_filename[sampid])
                         for sampid in sample_ids])

    bam_index_filename = dict([(sampid, bam_index_filename[sampid])
                         for sampid in sample_ids])

    gc_metrics = dict([(sampid, gc_metrics[sampid])
                        for sampid in sample_ids])

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )


    workflow.transform(
        name='postprocess_bam',
        ctx={'mem': config["memory"]['med']},
        axes=('sample_id',),
        func=tasks.postprocess_bam,
        args=(
            mgd.InputFile('merged_realign.bam', 'sample_id', fnames=bam),
            mgd.OutputFile('sorted_markdups', 'sample_id', fnames=bam_filename),
            mgd.OutputFile('sorted_markdups_index', 'sample_id', fnames=bam_index_filename),
            mgd.TempSpace('tempdir', 'sample_id'),
            config,
            mgd.OutputFile(markdups_metrics_filename, 'sample_id'),
            mgd.OutputFile(flagstat_metrics_filename, 'sample_id'),
        ),
    )
    
    workflow.transform(
        name='bam_collect_wgs_metrics',
        ctx={'mem': config["memory"]['med']},
        func=tasks.bam_collect_wgs_metrics,
        axes=('sample_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'sample_id', fnames=bam_filename),
            ref_genome,
            mgd.OutputFile(wgs_metrics_filename, 'sample_id'),
            config,
            mgd.TempSpace('wgs_tempdir', 'sample_id'),
        ),
    )
    
    workflow.transform(
        name='bam_collect_gc_metrics',
        ctx={'mem': config["memory"]['med']},
        func=tasks.bam_collect_gc_metrics,
        axes=('sample_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'sample_id', fnames=bam_filename),
            ref_genome,
            mgd.OutputFile(gc_metrics_filename, 'sample_id'),
            mgd.OutputFile(gc_summary_filename, 'sample_id'),
            mgd.OutputFile(gc_chart_filename, 'sample_id'),
            mgd.TempSpace('gc_tempdir', 'sample_id'),
        ),
    )
    
    workflow.transform(
        name='bam_collect_insert_metrics',
        ctx={'mem': config["memory"]['med']},
        func=tasks.bam_collect_insert_metrics,
        axes=('sample_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'sample_id', fnames=bam_filename),
            mgd.InputFile(flagstat_metrics_filename, 'sample_id'),
            mgd.OutputFile(insert_metrics_filename, 'sample_id'),
            mgd.OutputFile(insert_histogram_filename, 'sample_id'),
            mgd.TempSpace('insert_tempdir', 'sample_id'),
        ),
    )

    workflow.transform(
        name="collect_gc_metrics",
        func=tasks.collect_gc,
        ctx={'mem': config["memory"]['med']},
        axes = ('sample_id',),
        args = (
            mgd.InputFile(gc_metrics_filename, 'sample_id'),
            mgd.OutputFile('gc_metrics.csv', 'sample_id', fnames = gc_metrics),
            mgd.InputInstance('sample_id')
        )
    )
        
    workflow.transform(
        name='collect_metrics',
        ctx={'mem': config["memory"]['low']},
        func=tasks.collect_metrics,
        args=(
            mgd.InputFile(flagstat_metrics_filename, 'sample_id', axes_origin=[]),
            mgd.InputFile(markdups_metrics_filename, 'sample_id', axes_origin=[]),
            mgd.InputFile(insert_metrics_filename, 'sample_id', axes_origin=[]),
            mgd.InputFile(wgs_metrics_filename, 'sample_id', axes_origin=[]),
            mgd.TempOutputFile('metrics_summary.csv', 'sample_id', axes_origin=[]),
            mgd.OutputFile(alignment_metrics),
        ),
    )
    
    return workflow
