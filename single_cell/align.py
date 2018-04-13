'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
from workflows import alignment
from single_cell.utils import helpers


def align_workflow(workflow, args):

    args = helpers.generate_configs_in_temp(args)

    config = helpers.load_config(args)

    fastq1_files, fastq2_files  = helpers.get_fastqs(args['input_yaml'])
    seqinfo = helpers.get_seqinfo(args['input_yaml'])
    sampleinfo = helpers.get_sample_info(args['input_yaml'])

    cellids = helpers.get_samples(args['input_yaml'])
    bam_files, bai_files  = helpers.get_bams(args['input_yaml'])

    alignment_metrics = os.path.join(args["out_dir"], "metrics",'alignment_metrics.txt')
    gc_metrics = os.path.join(args["out_dir"], "metrics",'gc_metrics.txt')

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id', 'lane'),
        value=fastq1_files.keys(),
    )

    workflow.subworkflow(
        name='alignment_workflow',
        func=alignment.create_alignment_workflow,
        args=(
            mgd.InputFile('fastq_1', 'cell_id', 'lane', fnames=fastq1_files, axes_origin=[]),
            mgd.InputFile('fastq_2', 'cell_id', 'lane', fnames=fastq2_files, axes_origin=[]),
            mgd.OutputFile('bam_markdups', 'cell_id', fnames = bam_files, axes_origin=[]),
            mgd.OutputFile('bai_markdups', 'cell_id', fnames = bai_files, axes_origin=[]),
            mgd.OutputFile(alignment_metrics),
            mgd.OutputFile(gc_metrics),
            config['ref_genome'],
            config,
            args,
            seqinfo,
            sampleinfo,
            cellids,
        ),
    )

    return workflow
