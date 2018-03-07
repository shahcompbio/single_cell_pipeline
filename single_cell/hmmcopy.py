'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
from workflows import hmmcopy 
from utils import helpers

def hmmcopy_workflow(workflow, args):

    config = helpers.load_config(args)
    sampleids = helpers.get_samples(args['input_yaml'])
    bam_files, bai_files  = helpers.get_bams(args['input_yaml'])

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sampleids,
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
            args=(mgd.InputFile('bam_markdups', 'sample_id', fnames=bam_files),
                mgd.InputFile('bam_markdups_index', 'sample_id', fnames=bai_files),
                mgd.OutputFile(reads_filename),
                mgd.OutputFile(segs_filename),
                mgd.OutputFile(metrics_filename),
                mgd.OutputFile(params_filename),
                sampleids,
                config,
                args,
                params,
                results_dir
            ),
        )

    return workflow