'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import re

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers

from workflows import align, hmmcopy, qc_annotation


def qc_workflow(args):
    config = helpers.load_config(args)

    sampleinfo = helpers.get_sample_info(args['input_yaml'])
    cellids = helpers.get_samples(args['input_yaml'])
    bam_files, _ = helpers.get_bams(args['input_yaml'])

    lib = args["library_id"]

    workflow = pypeliner.workflow.Workflow()

    annotation_only = args['annotation_only']

    alignment_dir = args["alignment_output"]
    hmmcopy_dir = args["hmmcopy_output"]
    annotation_dir = args["annotation_output"]

    if alignment_dir and not annotation_only:
        alignment_files = get_output_files(alignment_dir, 'alignment', lib)

        fastq1_files, fastq2_files = helpers.get_fastqs(args['input_yaml'])
        triminfo = helpers.get_trim_info(args['input_yaml'])
        centerinfo = helpers.get_center_info(args['input_yaml'])

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
            kwargs={'realign': args['realign']}
        )

    if hmmcopy_dir and not annotation_only:
        hmmcopy_files = get_output_files(hmmcopy_dir, 'hmmcopy', lib)

        if not alignment_dir:
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

    if annotation_dir:
        annotation_files = get_output_files(annotation_dir, 'annotation', lib)
        if not hmmcopy_dir or not alignment_dir:
            raise Exception('--hmmcopy_output and --alignment_output are required to run annotation')

        alignment_files = get_output_files(alignment_dir, 'alignment', lib)
        hmmcopy_files = get_output_files(hmmcopy_dir, 'hmmcopy', lib)

        workflow.subworkflow(
            name='annotation_workflow',
            ctx={'docker_image': config['annotation']['docker']['single_cell_pipeline']},
            func=qc_annotation.create_qc_annotation_workflow,
            args=(
                mgd.InputFile(hmmcopy_files['metrics_csvs']),
                mgd.InputFile(hmmcopy_files['reads_csvs']),
                mgd.InputFile(alignment_files['alignment_metrics_csv']),
                mgd.InputFile(alignment_files['gc_metrics_csv']),
                mgd.InputFile(hmmcopy_files['segs_pdf']),
                mgd.OutputFile(annotation_files['merged_metrics_csvs']),
                mgd.OutputFile(annotation_files['qc_report']),
                mgd.OutputFile(annotation_files['corrupt_tree_newick']),
                mgd.OutputFile(annotation_files['consensus_tree_newick']),
                mgd.OutputFile(annotation_files['phylo_csv']),
                mgd.OutputFile(annotation_files['loci_rank_trees']),
                mgd.OutputFile(annotation_files['filtered_data']),
                mgd.OutputFile(annotation_files['corrupt_tree_pdf']),
                mgd.OutputFile(annotation_files['segs_pass']),
                mgd.OutputFile(annotation_files['segs_fail']),
                mgd.OutputFile(annotation_files['corrupt_heatmap_pdf']),
                mgd.OutputFile(annotation_files['heatmap_filt_pdf']),
                config['annotation'],
                lib,
            ),
            kwargs={'no_corrupt_tree': args['no_corrupt_tree']}
        )

    return workflow


def get_output_files(outdir, pipeline_type, lib):
    if pipeline_type == 'annotation':
        data = {
            'merged_metrics_csvs': os.path.join(outdir, '{0}_metrics.csv.gz'.format(lib)),
            'qc_report': os.path.join(outdir, '{0}_QC_report.html'.format(lib)),
            'corrupt_tree_newick': os.path.join(outdir, '{0}_corrupt_tree.newick'.format(lib)),
            'consensus_tree_newick': os.path.join(outdir, '{0}_corrupt_tree_consensus.newick'.format(lib)),
            'phylo_csv': os.path.join(outdir, '{0}_phylo.csv'.format(lib)),
            'loci_rank_trees': os.path.join(outdir, '{0}_rank_loci_trees.csv'.format(lib)),
            'filtered_data': os.path.join(outdir, '{0}_filtered_data.csv'.format(lib)),
            'corrupt_tree_pdf': os.path.join(outdir, '{0}_corrupt_tree.pdf'.format(lib)),
            'segs_pass': os.path.join(outdir, '{0}_segs_pass.tar.gz'.format(lib)),
            'segs_fail': os.path.join(outdir, '{0}_segs_fail.tar.gz'.format(lib)),
            'corrupt_heatmap_pdf': os.path.join(outdir, '{}_heatmap_corrupt_tree.pdf'.format(lib)),
            'heatmap_filt_pdf': os.path.join(outdir, '{}_heatmap_by_ec_filtered.pdf'.format(lib)),
        }

    elif pipeline_type == 'hmmcopy':
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

    elif pipeline_type == 'alignment':
        data = {
            'alignment_metrics_csv': os.path.join(outdir, '{}_alignment_metrics.csv.gz'.format(lib)),
            'gc_metrics_csv': os.path.join(outdir, '{}_gc_metrics.csv.gz'.format(lib)),
            'fastqc_metrics_csv': os.path.join(outdir, '{}_detailed_fastqscreen_metrics.csv.gz'.format(lib)),
            'plot_metrics_output': os.path.join(outdir, '{}_alignment_metrics.pdf'.format(lib)),
            'alignment_metrics_tar': os.path.join(outdir, '{}_alignment_metrics.tar.gz'.format(lib)),
        }
    else:
        raise Exception("Unknown pipeline type: {}".format(pipeline_type))

    return data


def generate_meta_files(args):
    alignment_dir = args["alignment_output"]
    hmmcopy_dir = args["hmmcopy_output"]
    annotation_dir = args["annotation_output"]
    annotation_only = args['annotation_only']
    lib = args["library_id"]

    cellids = helpers.get_samples(args['input_yaml'])

    samples = [re.split('[_-]', cell)[0] for cell in cellids]
    samples = sorted(set(samples))
    metadata = {
        'library_id': lib,
        'sample_ids': samples,
    }

    if annotation_dir:
        if not hmmcopy_dir or not alignment_dir:
            raise Exception('--hmmcopy_output and --alignment_output are required to run annotation')

    if alignment_dir and not annotation_only:
        alignment_files = get_output_files(alignment_dir, 'alignment', lib)
        metadata['type'] = 'align'
        helpers.generate_and_upload_metadata(
            args,
            alignment_dir,
            alignment_files.values(),
            metadata,
        )

    if hmmcopy_dir and not annotation_only:
        hmmcopy_files = get_output_files(hmmcopy_dir, 'hmmcopy', lib)
        metadata['type'] = 'hmmcopy'
        helpers.generate_and_upload_metadata(
            args,
            hmmcopy_dir,
            hmmcopy_files.values(),
            metadata,
        )

    if annotation_dir:
        annotation_files = get_output_files(annotation_dir, 'annotation', lib)
        metadata['type'] = 'annotation'
        helpers.generate_and_upload_metadata(
            args,
            annotation_dir,
            annotation_files.values(),
            metadata,
        )


def qc_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = qc_workflow(args)

    pyp.run(workflow)

    generate_meta_files(args)
