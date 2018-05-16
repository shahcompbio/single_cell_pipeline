'''
Created on Feb 20, 2018

@author: dgrewal
'''


import os
import pypeliner.managed as mgd
from workflows import aneufinder 
from single_cell.utils import helpers


def aneufinder_workflow(workflow, args):

    config = helpers.load_config(args)
    cellids = helpers.get_samples(args['input_yaml'])
    bam_files, _  = helpers.get_bams(args['input_yaml'])

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cellids,
    )

    output = os.path.join(args['out_dir'], 'results', "aneufinder")

    helpers.makedirs(output)

    results_filename = os.path.join(output, '{}_results.h5'.format(args['library_id']))
    workflow.subworkflow(
        name='aneufinder_workflow',
        func=aneufinder.create_aneufinder_workflow,
        args=(
            mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_files),
            cellids,
            config,
            output,
            mgd.OutputFile(results_filename),
            args['library_id']
        ),
    )

    return workflow
