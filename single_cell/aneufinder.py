'''
Created on Feb 20, 2018

@author: dgrewal
'''


import os
import pypeliner
import pypeliner.managed as mgd
from workflows import aneufinder 
from single_cell.utils import helpers


def aneufinder_workflow(workflow, args):

    config = helpers.load_config(args)
    sampleids = helpers.get_samples(args['input_yaml'])
    bam_files, _  = helpers.get_bams(args['input_yaml'])

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sampleids,
    )

    output = os.path.join(args['out_dir'], 'results', "aneufinder")

    helpers.makedirs(output)

    segs_filename = os.path.join(output, '{}_segments.csv'.format(args['library_id']))
    reads_filename = os.path.join(output, '{}_reads.csv'.format(args['library_id']))
    workflow.subworkflow(
        name='aneufinder_workflow',
        func=aneufinder.create_aneufinder_workflow,
        args=(
            mgd.InputFile('bam_markdups', 'sample_id', fnames=bam_files),
            sampleids,
            config,
            output,
            mgd.OutputFile(segs_filename),
            mgd.OutputFile(reads_filename),
            args['library_id']
        ),
    )

    return workflow
