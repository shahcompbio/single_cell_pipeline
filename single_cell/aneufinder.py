'''
Created on Feb 20, 2018

@author: dgrewal
'''


import os
import pypeliner.managed as mgd
from workflows import aneufinder 
from single_cell.utils import helpers
import single_cell


def aneufinder_workflow(workflow, args):

    config = helpers.load_config(args)
    config = config['aneufinder']

    baseimage = config['docker']['single_cell_pipeline']

    cellids = helpers.get_samples(args['input_yaml'])
    bam_files, _ = helpers.get_bams(args['input_yaml'])

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cellids,
    )

    info_file = os.path.join(args["out_dir"],'results', 'aneufinder', "info.yaml")

    output = os.path.join(args['out_dir'], 'results', "aneufinder")

    aneufinder_pdf_file = os.path.join(
        output, 'plots', '{}_reads.pdf'.format(args['library_id']))

    helpers.makedirs(output)

    results_filename = os.path.join(output, '{}_results.h5'.format(args['library_id']))
    workflow.subworkflow(
        name='aneufinder_workflow',
        func=aneufinder.create_aneufinder_workflow,
        args=(
            mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_files),
            cellids,
            config,
            mgd.OutputFile(results_filename),
            mgd.OutputFile(aneufinder_pdf_file),
        ),
    )

    results = {
        'aneufinder_plot': helpers.format_file_yaml(aneufinder_pdf_file),
        'aneufinder_data':helpers.format_file_yaml(results_filename),
    }

    input_datasets = {k: helpers.format_file_yaml(v) for k, v in bam_files.iteritems()}

    metadata = {
        'aneufinder':{
            'reads_table': '/aneufinder/reads',
            'segments_table': '/aneufinder/segments/',
            'chromosomes': config['chromosomes'],
            'ref_genome': config['ref_genome'],
            'version': single_cell.__version__,
            'results': results,
            'containers': config['docker'],
            'input_datasets': input_datasets,
            'output_datasets': None
        }
    }

    workflow.transform(
        name='generate_meta_yaml',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.utils.helpers.write_to_yaml",
        args=(
            mgd.OutputFile(info_file),
            metadata
        )
    )

    return workflow
