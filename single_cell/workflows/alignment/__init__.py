'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import tasks
import pypeliner
import pypeliner.managed as mgd


def create_alignment_workflow(
        fastq_1_filename,
        fastq_2_filename,
        bam_filename,
        ref_genome,
        lane_id,
        sample_id,
        config,
        args,
        seqinfo):

    out_dir = args['out_dir']
    fastqc_dir = os.path.join(out_dir, 'fastqc', lane_id)

    metrics_dir = os.path.join(args['out_dir'], 'metrics_per_lane',
                               lane_id, 'flagstat')
    flagstat_metrics_filename = os.path.join(metrics_dir,'{}.txt'.format(sample_id))

    read_group = tasks.get_readgroup(lane_id, sample_id, args, config, seqinfo)

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
            name='trimfastqs',
            ctx={'mem': config['low_mem']},
            func=tasks.trim_fastqs,
            args=(
                  mgd.InputFile(fastq_1_filename),
                  mgd.InputFile(fastq_2_filename),
                  mgd.TempOutputFile('trim_r1.fastq.gz'),
                  mgd.TempOutputFile('trim_r2.fastq.gz'),
                  fastqc_dir,
                  sample_id,
                  mgd.TempSpace('fastqc_1_temp'),
                  seqinfo[sample_id],
                  config)
    )

    workflow.transform(
        name='align_reads',
        ctx={'mem': config['high_mem']},
        func=tasks.align_pe,
        args=(
            mgd.TempInputFile('trim_r1.fastq.gz'),
            mgd.TempInputFile('trim_r2.fastq.gz'),
            mgd.OutputFile(bam_filename),
            mgd.OutputFile(flagstat_metrics_filename),
            mgd.TempSpace('alignment_temp'),
            mgd.InputFile(ref_genome),
            config,
            read_group
        )
    )

    return workflow
