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


def create_wgs_workflow(
    bam,
    bam_filename,
    bam_index_filename,
    ref_genome,
    sample_ids,
    config,
    out_dir):
 
    metrics_dir = os.path.join(out_dir, 'pseudo_wgs')
    markdups_metrics_filename = os.path.join(metrics_dir, 'markdups_metrics', 'markdups_metrics.txt')


    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )

    workflow.transform(
        name='merge_bams',
        ctx={'mem': config["memory"]['high'], 'pool_id': config['pools']['highmem']},
        func=tasks.merge_bams,
        args=(
            mgd.InputFile('bam', 'sample_id', fnames=bam),
            mgd.TempOutputFile('merged.bam'),
            mgd.TempSpace("temp_wgs_merge_bams")
        ),
    )

    workflow.transform(
        name='bam_sort',
        ctx={'mem': config["memory"]['high'], 'pool_id': config['pools']['highmem']},
        func=tasks.bam_sort,
        args=(
            mgd.TempInputFile('merged.bam'),
            mgd.TempOutputFile('sorted.bam'),
            mgd.TempSpace("temp_wgs_sort_bam")
        ),
    )

    workflow.transform(
        name='bam_markdups',
        ctx={'mem': config["memory"]['high'], 'pool_id': config['pools']['highmem']},
        func=tasks.bam_markdups,
        args=(
            mgd.TempInputFile('sorted.bam'),
            mgd.OutputFile(bam_filename),
            mgd.OutputFile(markdups_metrics_filename),
            mgd.TempSpace("temp_wgs_mkdup_bam")
        ),
    )

    workflow.transform(
        name='bam_index',
        ctx={'mem': config["memory"]['low'], 'pool_id': config['pools']['standard']},
        func=tasks.index_bam,
        args=(
            mgd.InputFile(bam_filename),
            mgd.OutputFile(bam_index_filename),
        ),
    )

    return workflow
