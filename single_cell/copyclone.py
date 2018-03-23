'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import pypeliner.managed as mgd
from workflows import copyclone
from utils import helpers
from workflows import singlecell_summary

def copyclone_workflow(workflow, args):

    config = helpers.load_config(args)
    sampleids = helpers.get_samples(args['input_yaml'])
    bam_files, bai_files  = helpers.get_bams(args['input_yaml'])

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sampleids,
    )

    results_dir = os.path.join(args['out_dir'], 'results', "copyclone")
    segs_filename = os.path.join(results_dir, '{}_segments.csv'.format(args['library_id']))
    reads_filename = os.path.join(results_dir, '{}_reads.csv'.format(args['library_id']))
    metrics_filename = os.path.join(results_dir, '{}_metrics.csv'.format(args['library_id']))


    workflow.subworkflow(
        name='copyclone_workflow',
        func=copyclone.create_copyclone_workflow,
        args=(mgd.InputFile('bam_markdups', 'sample_id', fnames=bam_files),
            mgd.InputFile('bam_markdups_index', 'sample_id', fnames=bai_files),
            mgd.OutputFile(reads_filename),
            mgd.OutputFile(segs_filename),
            mgd.OutputFile(metrics_filename),
            sampleids,
            config,
            args,
            results_dir
        ),
    )

    #assume these exist from prev runs for now
    alignment_metrics = os.path.join(args["out_dir"], "metrics",'alignment_metrics.txt')
    gc_metrics = os.path.join(args["out_dir"], "metrics",'gc_metrics.txt')

    reads_pdf_output = os.path.join(results_dir, 'plots', '{}_reads.pdf'.format(args['library_id']))
    segs_pdf_output = os.path.join(results_dir, 'plots', '{}_segs.pdf'.format(args['library_id']))
    bias_pdf_output = os.path.join(results_dir, 'plots', '{}_bias.pdf'.format(args['library_id']))

    params = {"num_states":7}

    # merge all samples per lane together
    workflow.subworkflow(
        name='summary_workflow_copyclone',
        func=singlecell_summary.create_summary_workflow,
        args=(
            mgd.InputFile(alignment_metrics),
            mgd.InputFile(gc_metrics),
            mgd.InputFile(segs_filename),
            mgd.InputFile(reads_filename),
            mgd.InputFile(metrics_filename),
            None,
            mgd.OutputFile(reads_pdf_output),
            mgd.OutputFile(segs_pdf_output),
            mgd.OutputFile(bias_pdf_output),
            None,
            config,
            params,
            results_dir,
            args,
            sampleids
        ),
    )

    return workflow