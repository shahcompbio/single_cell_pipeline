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

    info_file = os.path.join(args["out_dir"],'results','ltm', "info.yaml")

    results = {
        'ltm_cn_matrix': helpers.format_file_yaml(cn_matrix),
        'ltm_gml': helpers.format_file_yaml(output_gml),
        'ltm_rooted_gml': helpers.format_file_yaml(output_rooted_gml),
        'ltm_cnv_annots_csv': helpers.format_file_yaml(cnv_annots_csv),
        'ltm_cnv_tree_edges_csv': helpers.format_file_yaml(cnv_tree_edges_csv),
        'ltm_cnv_data_csv': helpers.format_file_yaml(cnv_data_csv),
        'ltm_output_rmd': helpers.format_file_yaml(output_rmd)
    }

    input_datasets = {k: helpers.format_file_yaml(v) for k,v in bam_file.iteritems()}

    metadata = {
        'LTM':{
            'chromosomes': config['chromosomes'],
            'ref_genome': config['ref_genome'],
            'cell_filters': config["good_cells"],
            'version': single_cell.__version__,
            'results': results,
            'containers': config['containers'],
            'input_datasets': input_datasets,
            'output_datasets': None
        }
    }

    workflow.transform(
        name='generate_meta_yaml',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],),
        func="single_cell.utils.helpers.write_to_yaml",
        args=(
            mgd.OutputFile(info_file),
            metadata
        )
    )


    return workflow
