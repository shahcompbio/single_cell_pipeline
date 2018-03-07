'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
from workflows import singlecell_summary 
from single_cell.utils import helpers



def summary_workflow(workflow, args):

    config = helpers.load_config(args)
    sampleids = helpers.get_samples(args['input_yaml'])

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sampleids,
    )

    for name, params in config["hmmcopy_params"].iteritems():
        name = "hmmcopy_" + name
 
        #assume these exist from prev runs for now
        alignment_metrics = os.path.join(args["out_dir"], "metrics",'alignment_metrics.txt')
        gc_metrics = os.path.join(args["out_dir"], "metrics",'gc_metrics.txt')

        results_dir = os.path.join(args['out_dir'], 'results', name)
        segs_filename = os.path.join(results_dir, '{}_segments.csv'.format(args['library_id']))
        reads_filename = os.path.join(results_dir, '{}_reads.csv'.format(args['library_id']))
        metrics_filename = os.path.join(results_dir, '{}_metrics.csv'.format(args['library_id']))
        params_filename = os.path.join(results_dir, '{}_params.csv'.format(args['library_id']))

        reads_pdf_output = os.path.join(results_dir, 'plots', '{}_reads.pdf'.format(args['library_id']))
        segs_pdf_output = os.path.join(results_dir, 'plots', '{}_segs.pdf'.format(args['library_id']))
        bias_pdf_output = os.path.join(results_dir, 'plots', '{}_bias.pdf'.format(args['library_id']))
        params_pdf_output = os.path.join(results_dir, 'plots', '{}_params.pdf'.format(args['library_id']))

        # merge all samples per lane together
        workflow.subworkflow(
            name='summary_workflow_' + name,
            func=singlecell_summary.create_summary_workflow,
            args=(
                mgd.InputFile(alignment_metrics),
                mgd.InputFile(gc_metrics),
                mgd.InputFile(segs_filename),
                mgd.InputFile(reads_filename),
                mgd.InputFile(metrics_filename),
                mgd.InputFile(params_filename),
                mgd.InputFile(reads_pdf_output),
                mgd.InputFile(segs_pdf_output),
                mgd.InputFile(bias_pdf_output),
                mgd.InputFile(params_pdf_output),
                config,
                params,
                results_dir,
                args,
                sampleids
            ),
        )

    return workflow