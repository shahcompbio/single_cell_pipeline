'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd

import tasks


def create_alignment_workflow(
    fastq_1_filename,
    fastq_2_filename,
    bam_filename,
    ref_genome,
    library_id,
    lane_id,
    sample_id,
    config,
    args):

    read_group = tasks.get_readgroup(library_id, lane_id, sample_id, args, config)

    metrics_dir = os.path.join(args['out_dir'], 'metrics_per_lane', lane_id,)
    flagstat_metrics_filename = os.path.join(metrics_dir, '{}.flagstat_metrics.txt'.format(sample_id))

    workflow = pypeliner.workflow.Workflow()

    workflow.commandline(
        name='bwa_align',
        ctx={'mem': config['med_mem']},
        args=(
            mgd.InputFile(fastq_1_filename),
            mgd.InputFile(fastq_2_filename),
            mgd.TempOutputFile('aligned.bam'),
            mgd.TempSpace('alignment_temp'),
            mgd.InputFile(ref_genome),
            config,
            read_group
            )
        )

    workflow.transform(
        name='bam_sort',
        ctx={'mem': config['high_mem']},
        func=tasks.bam_sort,
        args=(
            mgd.TempInputFile('aligned.bam'),
            mgd.OutputFile(bam_filename),
            config
        ),
    )

    workflow.commandline(
        name='bam_flagstat',
        ctx={'mem': config['low_mem']},
        args=(
            config['samtools'], 'flagstat',
            mgd.InputFile(bam_filename),
            '>',
            mgd.OutputFile(flagstat_metrics_filename),
        ),
    )

    return workflow
