import os
import re

import pypeliner.managed as mgd
from single_cell.utils import inpututils
from single_cell.workflows import align

import pypeliner
import sys

def get_output_files(outdir, lib):
    data = {
        'alignment_metrics_csv': os.path.join(outdir, '{}_alignment_metrics.csv.gz'.format(lib)),
        'gc_metrics_csv': os.path.join(outdir, '{}_gc_metrics.csv.gz'.format(lib)),
        'fastqc_metrics_csv': os.path.join(outdir, '{}_detailed_fastqscreen_metrics.csv.gz'.format(lib)),
        'plot_metrics_output': os.path.join(outdir, '{}_alignment_metrics.pdf'.format(lib)),
        'alignment_metrics_tar': os.path.join(outdir, '{}_alignment_metrics.tar.gz'.format(lib)),
    }

    return data


def alignment_workflow(args):
    config = inpututils.load_config(args)
    config = config['alignment']

    lib = args["library_id"]
    alignment_dir = args["out_dir"]
    bams_dir = args["bams_dir"]

    sampleinfo = inpututils.get_sample_info(args['input_yaml'])
    laneinfo = inpututils.get_lane_info(args['input_yaml'])

    cellids = inpututils.get_samples(args['input_yaml'])
    fastq1_files, fastq2_files = inpututils.get_fastqs(args['input_yaml'])

    alignment_files = get_output_files(alignment_dir, lib)
    alignment_meta = os.path.join(alignment_dir, 'metadata.yaml')

    bam_files_template = os.path.join(bams_dir, '{cell_id}.bam')
    mt_bam_files_template = os.path.join(bams_dir, '{cell_id}_MT.bam')
    bams_meta = os.path.join(bams_dir, 'metadata.yaml')

    lanes = sorted(set([v[1] for v in fastq1_files.keys()]))
    cells = sorted(set([v[0] for v in fastq1_files.keys()]))

    input_yaml_blob = os.path.join(alignment_dir, 'input.yaml')

    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config['docker']['single_cell_pipeline']}
    )

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id', 'lane'),
        value=list(fastq1_files.keys()),
    )

    workflow.subworkflow(
        name='alignment_workflow',
        func=align.create_alignment_workflow,
        args=(
            mgd.InputFile('fastq_1', 'cell_id', 'lane', fnames=fastq1_files, axes_origin=[]),
            mgd.InputFile('fastq_2', 'cell_id', 'lane', fnames=fastq2_files, axes_origin=[]),
            mgd.OutputFile(
                'bam_markdups', 'cell_id', template=bam_files_template,
                axes_origin=[], extensions=['.bai']
            ),
            mgd.OutputFile(
                'mt_bam_markdups', 'cell_id', template=mt_bam_files_template,
                axes_origin=[], extensions=['.bai']
            ),
            mgd.OutputFile(alignment_files['alignment_metrics_csv']),
            mgd.OutputFile(alignment_files['gc_metrics_csv']),
            mgd.OutputFile(alignment_files['fastqc_metrics_csv']),
            mgd.OutputFile(alignment_files['plot_metrics_output']),
            config['ref_genome'],
            config,
            laneinfo,
            sampleinfo,
            cellids,
            mgd.OutputFile(alignment_files['alignment_metrics_tar']),
            lib,
        ),
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='single_cell.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            alignment_dir,
            list(alignment_files.values()),
            mgd.OutputFile(alignment_meta)
        ),
        kwargs={
            'input_yaml_data': inpututils.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {
                'library_id': lib,
                'cell_ids': cells,
                'lane_ids': lanes,
                'type': 'alignment'
            }
        }
    )

    workflow.transform(
        name='generate_meta_files_bams',
        func='single_cell.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            bams_dir,
            mgd.Template('aligned.bam', 'cell_id', template=bam_files_template),
            mgd.OutputFile(bams_meta)
        ),
        kwargs={
            'metadata': {
                'library_id': lib,
                'cell_ids': cells,
                'lane_ids': lanes,
                'type': 'cellbams'
            },
            'template': (mgd.InputChunks('cell_id'), bam_files_template, 'cell_id'),
        }
    )

    return workflow


def alignment_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = alignment_workflow(args)

    pyp.run(workflow)
