'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import sys
import re

import pypeliner.managed as mgd
from single_cell.utils import inpututils
from single_cell.workflows import hmmcopy

import pypeliner


def get_output_files(outdir, lib):
    data = {
        'reads_csvs': os.path.join(outdir, '{0}_reads.csv.gz'.format(lib)),
        'segs_csvs': os.path.join(outdir, '{0}_segments.csv.gz'.format(lib)),
        'params_csvs': os.path.join(outdir, '{0}_params.csv.gz'.format(lib)),
        'metrics_csvs': os.path.join(outdir, '{0}_hmmcopy_metrics.csv.gz'.format(lib)),
        'hmmcopy_data_tar': os.path.join(outdir, '{0}_hmmcopy_data.tar.gz'.format(lib)),
        'igv_csvs': os.path.join(outdir, '{0}_igv_segments.seg'.format(lib)),
        'segs_pdf': os.path.join(outdir, '{}_segs.tar.gz'.format(lib)),
        'bias_pdf': os.path.join(outdir, '{}_bias.tar.gz'.format(lib)),
        'heatmap_pdf': os.path.join(outdir, '{}_heatmap_by_ec.pdf'.format(lib)),
        'metrics_pdf': os.path.join(outdir, '{}_hmmcopy_metrics.pdf'.format(lib)),
        'kernel_density_pdf': os.path.join(outdir, '{}_kernel_density.pdf'.format(lib)),
    }

    return data


def hmmcopy_workflow(args):
    config = inpututils.load_config(args)
    config = config['hmmcopy']

    sampleinfo = inpututils.get_sample_info(args['input_yaml'])
    cellids = inpututils.get_samples(args['input_yaml'])
    bam_files = inpututils.get_bams(args['input_yaml'])

    lib = args["library_id"]

    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config['docker']['single_cell_pipeline']},
    )

    hmmcopy_dir = args["out_dir"]

    hmmcopy_files = get_output_files(hmmcopy_dir, lib)
    hmmcopy_meta = os.path.join(hmmcopy_dir, 'metadata.yaml')
    input_yaml_blob = os.path.join(hmmcopy_dir, 'input.yaml')

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=list(bam_files.keys()),
    )

    workflow.subworkflow(
        name='hmmcopy_workflow',
        func=hmmcopy.create_hmmcopy_workflow,
        args=(
            mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_files, extensions=['.bai']),
            mgd.OutputFile(hmmcopy_files['reads_csvs']),
            mgd.OutputFile(hmmcopy_files['segs_csvs']),
            mgd.OutputFile(hmmcopy_files['metrics_csvs']),
            mgd.OutputFile(hmmcopy_files['params_csvs']),
            mgd.OutputFile(hmmcopy_files['igv_csvs']),
            mgd.OutputFile(hmmcopy_files['segs_pdf']),
            mgd.OutputFile(hmmcopy_files['bias_pdf']),
            mgd.OutputFile(hmmcopy_files['heatmap_pdf']),
            mgd.OutputFile(hmmcopy_files['metrics_pdf']),
            mgd.OutputFile(hmmcopy_files['kernel_density_pdf']),
            mgd.OutputFile(hmmcopy_files['hmmcopy_data_tar']),
            cellids,
            config,
            sampleinfo
        ),
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='single_cell.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            hmmcopy_dir,
            list(hmmcopy_files.values()),
            mgd.OutputFile(hmmcopy_meta)
        ),
        kwargs={
            'input_yaml_data': inpututils.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {
                'library_id': lib,
                'cell_ids': list(bam_files.keys()),
                'type': 'hmmcopy',
            }
        }
    )

    return workflow


def hmmcopy_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = hmmcopy_workflow(args)

    pyp.run(workflow)
