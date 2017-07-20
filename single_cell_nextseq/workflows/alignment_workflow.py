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
    lane_id,
    sample_id,
    config,
    args, desc):

    read_group = tasks.get_readgroup(desc, lane_id, config, sample_id)

    metrics_dir = os.path.join(args['out_dir'], 'lanes', lane_id, 'metrics')
    flagstat_metrics_filename = os.path.join(metrics_dir,'flagstat_metrics',  '{}.flagstat_metrics.txt'.format(sample_id))

    workflow = pypeliner.workflow.Workflow()

    workflow.commandline(
        name='aln_read_1',
        ctx={'mem': config['med_mem']},
        args=(
            config['bwa'],
            'aln',
            mgd.InputFile(ref_genome),
            mgd.InputFile(fastq_1_filename),
            '>',
            mgd.TempOutputFile('read_1.sai'),
        ),
    )

    workflow.commandline(
        name='aln_read_2',
        ctx={'mem': config['med_mem']},
        args=(
            config['bwa'], 'aln',
            mgd.InputFile(ref_genome),
            mgd.InputFile(fastq_2_filename),
            '>',
            mgd.TempOutputFile('read_2.sai'),
        ),
    )

    workflow.commandline(
        name='sampe',
        ctx={'mem': config['med_mem']},
        args=(
            config['bwa'], 'sampe',
            '-r', read_group,
            mgd.InputFile(ref_genome),
            mgd.TempInputFile('read_1.sai'),
            mgd.TempInputFile('read_2.sai'),
            mgd.InputFile(fastq_1_filename),
            mgd.InputFile(fastq_2_filename),
            '|',
            config['samtools'], 'view',
            '-bSh', '-',
            '>',
            mgd.TempOutputFile('aligned.bam'),
        ),
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
