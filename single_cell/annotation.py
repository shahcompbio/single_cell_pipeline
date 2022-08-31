'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import sys

import pypeliner.managed as mgd
from single_cell.utils import inpututils
from single_cell.workflows import qc_annotation

import pypeliner


def annotation_workflow(args):
    config = inpututils.load_config(args)

    annotation_infiles = inpututils.load_yaml(args['input_yaml'])

    lib = args["library_id"]

    workflow = pypeliner.workflow.Workflow()

    annotation_dir = args["output_prefix"]

    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')
    annotation_files = get_output_files(annotation_dir)
    annotation_meta = os.path.join(args['out_dir'], 'metadata.yaml')

    workflow.subworkflow(
        name='annotation_workflow',
        func=qc_annotation.create_qc_annotation_workflow,
        args=(
            mgd.InputFile(annotation_infiles['hmmcopy_metrics']),
            mgd.InputFile(annotation_infiles['hmmcopy_reads']),
            mgd.InputFile(annotation_infiles['alignment_metrics']),
            mgd.InputFile(annotation_infiles['gc_metrics']),
            mgd.InputFile(annotation_infiles['segs_pdf_tar']),
            mgd.OutputFile(annotation_files['merged_metrics_csvs']),
            mgd.OutputFile(annotation_files['qc_report']),
            mgd.OutputFile(annotation_files['segs_pass']),
            mgd.OutputFile(annotation_files['segs_fail']),
            mgd.OutputFile(annotation_files['heatmap_filt_pdf']),
            config['annotation'],
        )
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='single_cell.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args['out_dir'],
            list(annotation_files.values()),
            mgd.OutputFile(annotation_meta)
        ),
        kwargs={
            'input_yaml_data': inpututils.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {
                'library_id': lib,
                'type': 'annotation'
            }
        }
    )

    return workflow


def get_output_files(outdir):
    data = {
        'merged_metrics_csvs': outdir + 'metrics.csv.gz',
        'qc_report': outdir + 'QC_report.html',
        'segs_pass': outdir + 'segs_pass.tar.gz',
        'segs_fail': outdir + 'segs_fail.tar.gz',
        'heatmap_filt_pdf': outdir + 'heatmap_by_ec_filtered.pdf',
    }

    return data


def annotation_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = annotation_workflow(args)

    pyp.run(workflow)
