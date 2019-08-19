'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import re

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers

from workflows import align, hmmcopy, qc_annotation


def qc_workflow(args):
    config = helpers.load_config(args)

    cellids = helpers.get_samples(args['input_yaml'])
    bam_files, _ = helpers.get_bams(args['input_yaml'])

    lib = args["library_id"]

    workflow = pypeliner.workflow.Workflow()

    annotation_only = args['annotation_only']

    alignment_dir = args["alignment_output"]

    if alignment_dir and not annotation_only:
        gc_file = os.path.join(alignment_dir, '{}_gc_metrics.csv.gz'.format(lib))

        fastq1_files, fastq2_files = helpers.get_fastqs(args['input_yaml'])

        workflow.setobj(
            obj=mgd.OutputChunks('cell_id', 'lane'),
            value=list(fastq1_files.keys()),
        )

        workflow.subworkflow(
            name='alignment_workflow',
            ctx={'docker_image': config['alignment']['docker']['single_cell_pipeline']},
            func=align.create_alignment_workflow,
            args=(
                mgd.OutputFile('bam_markdups', 'cell_id', fnames=bam_files, axes_origin=[], extensions=['.bai']),
                mgd.OutputFile(gc_file),
                config['alignment'],
                config['alignment']['ref_genome'],
                cellids,
            ),
        )


    return workflow

def qc_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = qc_workflow(args)

    pyp.run(workflow)

