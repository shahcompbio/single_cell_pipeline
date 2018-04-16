'''
Created on Apr 13, 2018

@author: dgrewal
'''

import os
import pypeliner
import pypeliner.managed as mgd
from workflows import titan
from single_cell.utils import helpers


def copy_number_calling_workflow(workflow, args):

    args = helpers.generate_configs_in_temp(args)

    config = helpers.load_config(args)

    tumour_bam_files, tumour_bai_files  = helpers.get_bams(args['tumour_yaml'])

    normal_bam_files, normal_bai_files  = helpers.get_bams(args['normal_yaml'])

    tumour_cellids = helpers.get_samples(args['tumour_yaml'])

    normal_cellids = helpers.get_samples(args['normal_yaml'])

    copynumber_dir = os.path.join(args["out_dir"], "copynumber")

    out_file = os.path.join(copynumber_dir,  "output.csv")

    workflow.setobj(
        obj=mgd.OutputChunks('tumour_cell_id'),
        value=tumour_cellids,
    )

    workflow.setobj(
        obj=mgd.OutputChunks('normal_cell_id'),
        value=tumour_cellids,
    )

    workflow.subworkflow(
        name='titan_workflow',
        func=titan.create_titan_workflow,
        args=(
            mgd.InputFile('bam_markdups', 'tumour_cell_id', fnames=tumour_bam_files),
            mgd.InputFile('bam_markdups_index', 'tumour_cell_id', fnames=tumour_bai_files),
            mgd.InputFile('bam_markdups', 'normal_cell_id', fnames=normal_bam_files),
            mgd.InputFile('bam_markdups_index', 'normal_cell_id', fnames=normal_bai_files),
            config['ref_genome'],
            copynumber_dir,
            out_file,
            config,
            args,
            tumour_cellids,
            normal_cellids,
        ),
    )


    return workflow