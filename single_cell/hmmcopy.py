'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import re

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers
from single_cell.utils import inpututils
from single_cell.workflows import hmmcopy


def hmmcopy_workflow(args):
    config = inpututils.load_config(args)

    sampleinfo = inpututils.get_sample_info(args['input_yaml'])
    cellids = inpututils.get_samples(args['input_yaml'])
    bam_files, _ = inpututils.get_bams(args['input_yaml'])

    lib = args["library_id"]

    workflow = pypeliner.workflow.Workflow()

    hmmcopy_dir = args["out_dir"]

    hmmcopy_files = get_output_files(hmmcopy_dir, lib)

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=list(bam_files.keys()),
    )

    workflow.subworkflow(
        name='hmmcopy_workflow',
        ctx={'docker_image': config['hmmcopy']['docker']['single_cell_pipeline']},
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
            config['hmmcopy'],
            sampleinfo
        ),
    )

    return workflow


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


def generate_meta_files(args):
    hmmcopy_dir = args["out_dir"]
    lib = args["library_id"]

    cellids = inpututils.get_samples(args['input_yaml'])

    samples = [re.split('[_-]', cell)[0] for cell in cellids]
    samples = sorted(set(samples))
    metadata = {
        'library_id': lib,
        'sample_ids': samples,
    }

    hmmcopy_files = get_output_files(hmmcopy_dir, lib)
    metadata['type'] = 'hmmcopy'
    helpers.generate_and_upload_metadata(
        args,
        hmmcopy_dir,
        hmmcopy_files.values(),
        metadata,
        input_yaml=args['input_yaml']
    )


def hmmcopy_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = hmmcopy_workflow(args)

    pyp.run(workflow)

    generate_meta_files(args)
