'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import sys

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import inpututils
from single_cell.workflows import hmmcopy


def get_output_files(outdir):
    data = {
        'reads_csvs': outdir + '{0}_reads.csv.gz',
        'segs_csvs': outdir + '_segments.csv.gz',
        'params_csvs': outdir + '_params.csv.gz',
        'metrics_csvs': outdir + '_hmmcopy_metrics.csv.gz',
        'hmmcopy_data_tar': outdir + '_hmmcopy_data.tar.gz',
        'igv_csvs': outdir + '_igv_segments.seg',
        'segs_pdf': outdir + '_segs.tar.gz',
        'bias_pdf': outdir + '_bias.tar.gz',
        'heatmap_pdf': outdir + '_heatmap_by_ec.pdf',
        'metrics_pdf': outdir + '_hmmcopy_metrics.pdf',
        'kernel_density_pdf': outdir + '_kernel_density.pdf',
    }

    return data


def hmmcopy_workflow(args):
    config = inpututils.load_config(args)
    config = config['hmmcopy']

    sampleinfo = inpututils.get_sample_info(args['input_yaml'])
    cellids = inpututils.get_samples(args['input_yaml'])
    bam_files = inpututils.get_bams(args['input_yaml'])

    lib = args["library_id"]

    workflow = pypeliner.workflow.Workflow()

    hmmcopy_prefix = args["output_prefix"]

    hmmcopy_files = get_output_files(hmmcopy_prefix)
    hmmcopy_meta = os.path.join(hmmcopy_prefix, 'metadata.yaml')
    input_yaml_blob = os.path.join(hmmcopy_prefix, 'input.yaml')

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
            hmmcopy_prefix,
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
