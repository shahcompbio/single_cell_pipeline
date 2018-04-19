'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
from workflows import hmmcopy 
from utils import helpers
from workflows import singlecell_summary 



def hmmcopy_workflow(workflow, args):

    config = helpers.load_config(args)
    cellids = helpers.get_samples(args['input_yaml'])
    bam_files, bai_files  = helpers.get_bams(args['input_yaml'])

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cellids,
    )

    for name, params in config["hmmcopy_params"].iteritems():
        name = "hmmcopy_" + name
        results_dir = os.path.join(args['out_dir'], 'results', name)
        segs_filename = os.path.join(results_dir, '{}_segments.csv'.format(args['library_id']))
        reads_filename = os.path.join(results_dir, '{}_reads.csv'.format(args['library_id']))
        metrics_filename = os.path.join(results_dir, '{}_metrics.csv'.format(args['library_id']))
        params_filename = os.path.join(results_dir, '{}_params.csv'.format(args['library_id']))

        workflow.subworkflow(
            name='hmmcopy_workflow_' + name,
            func=hmmcopy.create_hmmcopy_workflow,
            args=(mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_files),
                mgd.InputFile('bam_markdups_index', 'cell_id', fnames=bai_files),
                mgd.OutputFile(reads_filename),
                mgd.OutputFile(segs_filename),
                mgd.OutputFile(metrics_filename),
                mgd.OutputFile(params_filename),
                cellids,
                config,
                args,
                params,
                results_dir
            ),
        )

        #assume these exist from prev runs for now
        reads_pdf_output = os.path.join(results_dir, 'plots', '{}_reads.pdf'.format(args['library_id']))
        segs_pdf_output = os.path.join(results_dir, 'plots', '{}_segs.pdf'.format(args['library_id']))
        bias_pdf_output = os.path.join(results_dir, 'plots', '{}_bias.pdf'.format(args['library_id']))
        params_pdf_output = os.path.join(results_dir, 'plots', '{}_params.pdf'.format(args['library_id']))

        # merge all samples per lane together
        workflow.subworkflow(
            name='summary_workflow_' + name,
            func=singlecell_summary.create_summary_workflow,
            args=(
                mgd.InputFile(segs_filename),
                mgd.InputFile(reads_filename),
                mgd.InputFile(metrics_filename),
                mgd.InputFile(params_filename),
                mgd.OutputFile(reads_pdf_output),
                mgd.OutputFile(segs_pdf_output),
                mgd.OutputFile(bias_pdf_output),
                mgd.OutputFile(params_pdf_output),
                config,
                params,
                results_dir,
                args,
                cellids
            ),
        )



    return workflow