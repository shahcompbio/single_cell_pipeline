import os
import re

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers
from single_cell.utils import inpututils
from single_cell.workflows import align


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

    sampleinfo = inpututils.get_sample_info(args['input_yaml'])
    cellids = inpututils.get_samples(args['input_yaml'])
    bam_files, _ = inpututils.get_bams(args['input_yaml'])

    lib = args["library_id"]

    workflow = pypeliner.workflow.Workflow()

    alignment_dir = args["out_dir"]
    alignment_files = get_output_files(alignment_dir, lib)

    fastq1_files, fastq2_files = inpututils.get_fastqs(args['input_yaml'])
    triminfo = inpututils.get_trim_info(args['input_yaml'])
    centerinfo = inpututils.get_center_info(args['input_yaml'])

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id', 'lane'),
        value=list(fastq1_files.keys()),
    )

    workflow.subworkflow(
        name='alignment_workflow',
        ctx={'docker_image': config['alignment']['docker']['single_cell_pipeline']},
        func=align.create_alignment_workflow,
        args=(
            mgd.InputFile('fastq_1', 'cell_id', 'lane', fnames=fastq1_files, axes_origin=[]),
            mgd.InputFile('fastq_2', 'cell_id', 'lane', fnames=fastq2_files, axes_origin=[]),
            mgd.OutputFile('bam_markdups', 'cell_id', fnames=bam_files, axes_origin=[], extensions=['.bai']),
            mgd.OutputFile(alignment_files['alignment_metrics_csv']),
            mgd.OutputFile(alignment_files['gc_metrics_csv']),
            mgd.OutputFile(alignment_files['fastqc_metrics_csv']),
            mgd.OutputFile(alignment_files['plot_metrics_output']),
            config['alignment']['ref_genome'],
            config['alignment'],
            triminfo,
            centerinfo,
            sampleinfo,
            cellids,
            mgd.OutputFile(alignment_files['alignment_metrics_tar']),
            lib,
        ),
    )

    return workflow


def generate_meta_files(args):
    alignment_dir = args["out_dir"]
    lib = args["library_id"]

    cellids = inpututils.get_samples(args['input_yaml'])

    samples = [re.split('[_-]', cell)[0] for cell in cellids]
    samples = sorted(set(samples))
    metadata = {
        'library_id': lib,
        'sample_ids': samples,
    }

    alignment_files = get_output_files(alignment_dir, lib)
    metadata['type'] = 'align'
    helpers.generate_and_upload_metadata(
        args,
        alignment_dir,
        alignment_files.values(),
        metadata,
        input_yaml=args['input_yaml']
    )


def alignment_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = alignment_workflow(args)

    pyp.run(workflow)

    generate_meta_files(args)
