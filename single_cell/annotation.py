
'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import re

import pypeliner.managed as mgd
from single_cell.utils import helpers
from single_cell.utils import inpututils
from single_cell.workflows import align, hmmcopy, qc_annotation

import pypeliner


def annotation_workflow(args):
    config = inpututils.load_config(args)

    bam_files, _ = inpututils.get_bams(args['input_yaml'])

    lib = args["library_id"]

    workflow = pypeliner.workflow.Workflow()

    annotation_dir = args["annotation_output"]

    annotation_files = get_output_files(annotation_dir, 'annotation', lib)

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


def get_output_files(outdir, lib):
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

    return data


def generate_meta_files(args):
    annotation_dir = args["out_dir"]
    lib = args["library_id"]

    cellids = inpututils.get_samples(args['input_yaml'])

    samples = [re.split('[_-]', cell)[0] for cell in cellids]
    samples = sorted(set(samples))
    metadata = {
        'library_id': lib,
        'sample_ids': samples,
    }

    annotation_files = get_output_files(annotation_dir, lib)
    metadata['type'] = 'annotation'
    helpers.generate_and_upload_metadata(
        args,
        annotation_dir,
        annotation_files.values(),
        metadata,
        input_yaml=args['input_yaml']
    )


def annotation_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = annotation_workflow(args)

    pyp.run(workflow)

    generate_meta_files(args)
