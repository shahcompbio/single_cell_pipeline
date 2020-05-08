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

    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config['annotation']['docker']['single_cell_pipeline']},
    )

    annotation_dir = args["out_dir"]

    input_yaml_blob = os.path.join(annotation_dir, 'input.yaml')
    annotation_files = get_output_files(annotation_dir, lib)
    annotation_meta = os.path.join(annotation_dir, 'metadata.yaml')

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
            config['alignment'],
        )
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='single_cell.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            annotation_dir,
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


def get_output_files(outdir, lib):
    data = {
        'merged_metrics_csvs': os.path.join(outdir, '{0}_metrics.csv.gz'.format(lib)),
        'qc_report': os.path.join(outdir, '{0}_QC_report.html'.format(lib)),
        'segs_pass': os.path.join(outdir, '{0}_segs_pass.tar.gz'.format(lib)),
        'segs_fail': os.path.join(outdir, '{0}_segs_fail.tar.gz'.format(lib)),
        'heatmap_filt_pdf': os.path.join(outdir, '{}_heatmap_by_ec_filtered.pdf'.format(lib)),
    }

    return data


def annotation_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = annotation_workflow(args)

    pyp.run(workflow)
