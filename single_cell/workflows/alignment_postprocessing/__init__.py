'''
Created on Jul 11, 2017

@author: dgrewal
'''


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

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )


    workflow.transform(
        name='bam_sort',
        ctx={'mem': config['high_mem']},
        axes=('sample_id',),
        func=tasks.bam_sort,
        args=(
            mgd.InputFile('merged_realign.bam', 'sample_id', fnames=bam),
            mgd.TempOutputFile('sorted.bam', 'sample_id'),
            config
        ),
    )
   
    workflow.transform(
        name='bam_markdups',
        ctx={'mem': config['high_mem']},
        func=tasks.bam_markdups,
        args=(
            mgd.TempInputFile('sorted.bam', 'sample_id'),
            mgd.OutputFile('sorted_markdups', 'sample_id', fnames=bam_filename),
            mgd.OutputFile(markdups_metrics_filename, 'sample_id'),
            config
        ),
    )
   
    workflow.commandline(
        name='bam_index',
        ctx={'mem': config['low_mem']},
        axes=('sample_id',),
        args=(
            'samtools', 'index',
            mgd.InputFile('sorted_markdups', 'sample_id', fnames=bam_filename),
            mgd.OutputFile('sorted_markdups_index', 'sample_id', fnames=bam_index_filename),
        ),
    )
   
    workflow.commandline(
        name='bam_flagstat',
        ctx={'mem': config['low_mem']},
        axes=('sample_id',),
        args=(
            'samtools', 'flagstat',
            mgd.InputFile('sorted_markdups', 'sample_id', fnames=bam_filename),
            '>',
            mgd.OutputFile(flagstat_metrics_filename, 'sample_id'),
        ),
    )
   
    workflow.transform(
        name='bam_collect_wgs_metrics',
        ctx={'mem': config['high_mem']},
        func=tasks.bam_collect_wgs_metrics,
        axes=('sample_id',),
        args=(
            mgd.InputFile(bam_filename),
            mgd.InputFile(ref_genome),
            mgd.OutputFile(wgs_metrics_filename, 'sample_id'),
            config,
        ),
    )
   
    workflow.transform(
        name='bam_collect_gc_metrics',
        ctx={'mem': config['high_mem']},
        func=tasks.bam_collect_gc_metrics,
        axes=('sample_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'sample_id', fnames=bam_filename),
            mgd.InputFile(ref_genome),
            mgd.OutputFile(gc_metrics_filename, 'sample_id'),
            mgd.OutputFile(gc_summary_filename, 'sample_id'),
            mgd.OutputFile(gc_chart_filename, 'sample_id'),
            config
        ),
    )
   
    workflow.transform(
        name='bam_collect_insert_metrics',
        ctx={'mem': config['high_mem']},
        func=tasks.bam_collect_insert_metrics,
        axes=('sample_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'sample_id', fnames=bam_filename),
            mgd.InputFile(flagstat_metrics_filename, 'sample_id'),
            mgd.OutputFile(insert_metrics_filename, 'sample_id'),
            mgd.OutputFile(insert_histogram_filename, 'sample_id'),
            config
        ),
    )
       
    workflow.transform(
        name='collect_metrics',
        ctx={'mem': config['low_mem']},
        func=tasks.collect_metrics,
        axes=('sample_id',),
        args=(
            mgd.InputFile(flagstat_metrics_filename, 'sample_id'),
            mgd.InputFile(markdups_metrics_filename, 'sample_id'),
            mgd.InputFile(insert_metrics_filename, 'sample_id'),
            mgd.InputFile(wgs_metrics_filename, 'sample_id'),
            mgd.TempOutputFile('metrics_summary.csv', 'sample_id'),
            mgd.InputInstance('sample_id'),
        ),
    )
   
    workflow.transform(
        name='collect_gc_metrics',
        ctx={'mem': config['low_mem']},
        func = tasks.collect_gc_metrics,
        axes=('sample_id',),
        args=(
            mgd.InputFile(gc_metrics_filename, 'sample_id'),
            mgd.TempOutputFile('gc_matrix.csv', 'sample_id'),
            ',',
            'NORMALIZED_COVERAGE',
            mgd.InputInstance('sample_id'),
            'gcbias'
        ),
    )


    workflow.transform(
        name='merge_summary_metrics',
        ctx={'mem': config['low_mem']},
        func=tasks.concatenate_csv,
        args=(
            mgd.TempInputFile('metrics_summary.csv', 'sample_id'),
            mgd.OutputFile(alignment_metrics),
        ),
    )

    workflow.transform(
        name='merge_gc_metrics',
        ctx={'mem': config['low_mem']},
        func=tasks.merge_csv,
        args=(
            mgd.TempInputFile('gc_matrix.csv', 'sample_id'),
            mgd.OutputFile(gc_metrics),
            'outer',
            'gc'
        ),
    )

    return workflow
