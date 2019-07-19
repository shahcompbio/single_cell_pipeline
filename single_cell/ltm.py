'''
Created on July 30, 2018

@author: pwalters
'''

import os
import pypeliner
import pypeliner.managed as mgd

from workflows import ltm
from single_cell.utils import helpers
from single_cell.utils import ltmutils
import single_cell


def ltm_workflow(args):
    workflow = pypeliner.workflow.Workflow()

    config = helpers.load_config(args)

    hmmcopy, timepoints = ltmutils.read_input_file(args['input_csv'])

    cn_matrix = os.path.join(args['out_dir'], 'cn_matrix.csv')
    output_gml = os.path.join(args['out_dir'], 'tree.gml')
    output_rooted_gml = os.path.join(args['out_dir'], 'rooted_tree.gml')

    # Outputs required for visualization with cellscape
    cnv_annots_csv = os.path.join(args['out_dir'], 'cnv_annots.csv')
    cnv_tree_edges_csv = os.path.join(args['out_dir'], 'cnv_tree_edges.csv')
    cnv_data_csv = os.path.join(args['out_dir'], 'cnv_data.csv')
    output_rmd = os.path.join(args['out_dir'], 'cellscape.Rmd')
    root_id_file = os.path.join(args['out_dir'], 'root_id.txt')

    workflow.setobj(
        obj=mgd.OutputChunks('timepoint'),
        value=timepoints,
    )

    workflow.subworkflow(
        name='ltm_scale',
        func=ltm.create_ltm_workflow,
        args=(
            mgd.InputFile('hmmcopy.h5', 'timepoint', fnames=hmmcopy),
            mgd.OutputFile(cn_matrix),
            mgd.OutputFile(output_gml),
            mgd.OutputFile(output_rooted_gml),
            mgd.OutputFile(cnv_annots_csv),
            mgd.OutputFile(cnv_tree_edges_csv),
            mgd.OutputFile(cnv_data_csv),
            mgd.OutputFile(output_rmd),
            config,
            args['root_id'],
            mgd.OutputFile(root_id_file),
            args['number_of_jobs'],
            args['ploidy'],
        ),
    )

    return workflow

def ltm_pipeline(args):

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = ltm_workflow(args)

    pyp.run(workflow)