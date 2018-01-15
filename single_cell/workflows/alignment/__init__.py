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
        bam_filename,
        ref_genome,
        sample_lanes,
        config,
        args,
        seqinfo):

    out_dir = args['out_dir']

    fastqc_reports = os.path.join(out_dir, 'fastqc', "{lane}_{sample_id}_reports.tar.gz")


    metrics_dir = os.path.join(args['out_dir'], 'metrics_per_lane',
                               '{lane}', 'flagstat')
    flagstat_metrics_filename = os.path.join(metrics_dir,'{sample_id}.txt')

    bam_filename = dict([(key, bam_filename(key[0], key[1]))
                        for key in sample_lanes])


    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'lane'),
        value=fastq_1_filename.keys(),
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('seqinfo', 'sample_id', axes_origin=[]),
        value=seqinfo)

    workflow.transform(
        name='align_reads',
        ctx={'mem': config['med_mem']},
        axes = ('sample_id', 'lane',),
        func=tasks.align_pe,
        args=(
            mgd.InputFile('fastq_1', 'sample_id', 'lane', fnames=fastq_1_filename),
            mgd.OutputFile('aligned_per_cell_per_lane.sorted.bam', 'sample_id', 'lane', fnames=bam_filename),
            mgd.OutputFile(fastqc_reports, 'sample_id', 'lane'),
            mgd.OutputFile(flagstat_metrics_filename, 'sample_id', 'lane'),
            mgd.TempSpace('alignment_temp', 'sample_id', 'lane'),
            ref_genome,
            mgd.TempInputObj('seqinfo', 'sample_id'),
            mgd.InputInstance('sample_id'),
            mgd.InputInstance('lane'),
            args['library_id']
        )
    )

    return workflow
