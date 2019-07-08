'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers

from workflows import align, hmmcopy, qc_annotation


def qc_workflow(args):
    run_alignment = args['alignment']
    run_hmmcopy = args['hmmcopy']
    run_annotation = args['annotation']

    if not any((run_alignment, run_hmmcopy, run_annotation)):
        run_alignment = True
        run_hmmcopy = True
        run_annotation = True

    config = helpers.load_config(args)

    sampleinfo = helpers.get_sample_info(args['input_yaml'])
    cellids = helpers.get_samples(args['input_yaml'])
    bam_files, _ = helpers.get_bams(args['input_yaml'])

    lib = args["library_id"]

    outdir = os.path.join(args["out_dir"], "results", "QC")

    alignment_dir = os.path.join(outdir, 'alignment')
    alignment_metrics_csv = os.path.join(alignment_dir, '{}_alignment_metrics.csv.gz'.format(lib))
    gc_metrics_csv = os.path.join(alignment_dir, '{}_gc_metrics.csv.gz'.format(lib))
    fastqc_metrics_csv = os.path.join(alignment_dir, '{}_detailed_fastqscreen_metrics.csv.gz'.format(lib))
    plot_metrics_output = os.path.join(alignment_dir, '{}_plot_metrics.pdf'.format(lib))
    alignment_metrics_tar = os.path.join(alignment_dir, '{}_alignment_metrics.tar.gz'.format(lib))

    hmmcopy_dir = os.path.join(outdir, 'hmmcopy_autoploidy')
    reads_csvs = os.path.join(hmmcopy_dir, '{0}_reads.csv.gz'.format(lib))
    segs_csvs = os.path.join(hmmcopy_dir, '{0}_segments.csv.gz'.format(lib))
    params_csvs = os.path.join(hmmcopy_dir, '{0}_params.csv.gz'.format(lib))
    metrics_csvs = os.path.join(hmmcopy_dir, '{0}_metrics.csv.gz'.format(lib))
    hmmcopy_data_tar = os.path.join(hmmcopy_dir, '{0}_data.tar.gz'.format(lib))
    igv_csvs = os.path.join(hmmcopy_dir, '{0}_igv_segments.seg'.format(lib))
    segs_pdf = os.path.join(hmmcopy_dir, '{}_segs.tar.gz'.format(lib))
    bias_pdf = os.path.join(hmmcopy_dir, '{}_bias.tar.gz'.format(lib))
    heatmap_filt_pdf = os.path.join(hmmcopy_dir, '{}_heatmap_by_ec_filtered.pdf'.format(lib))
    heatmap_pdf = os.path.join(hmmcopy_dir, '{}_heatmap_by_ec.pdf'.format(lib))
    metrics_pdf = os.path.join(hmmcopy_dir, '{}_metrics.pdf'.format(lib))
    kernel_density_pdf = os.path.join(hmmcopy_dir, '{}_kernel_density.pdf'.format(lib))

    annotation_dir = os.path.join(outdir, 'annotation')
    merged_metrics_csvs = os.path.join(annotation_dir, '{0}_metrics.csv.gz'.format(lib))
    qc_report = os.path.join(annotation_dir, '{0}_QC_report.html'.format(lib))
    corrupt_tree_newick = os.path.join(annotation_dir, '{0}_corrupt_tree.newick'.format(lib))
    consensus_tree_newick = os.path.join(annotation_dir, '{0}_corrupt_tree_consensus.newick'.format(lib))
    phylo_csv = os.path.join(annotation_dir, '{0}_phylo.csv'.format(lib))
    loci_rank_trees = os.path.join(annotation_dir, '{0}_rank_loci_trees.csv'.format(lib))
    filtered_data = os.path.join(annotation_dir, '{0}_filtered_data.csv'.format(lib))
    corrupt_tree_pdf = os.path.join(annotation_dir, '{0}_corrupt_tree.pdf'.format(lib))
    segs_pass = os.path.join(annotation_dir, '{0}_segs_pass.tar.gz'.format(lib))
    segs_fail = os.path.join(annotation_dir, '{0}_segs_fail.tar.gz'.format(lib))
    corrupt_heatmap_pdf = os.path.join(annotation_dir, '{}_heatmap_corrupt_tree.pdf'.format(lib))

    workflow = pypeliner.workflow.Workflow()

    if run_alignment:
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
                mgd.OutputFile(alignment_metrics_csv),
                mgd.OutputFile(gc_metrics_csv),
                mgd.OutputFile(fastqc_metrics_csv),
                mgd.OutputFile(plot_metrics_output),
                config['alignment']['ref_genome'],
                config['alignment'],
                triminfo,
                centerinfo,
                sampleinfo,
                cellids,
                mgd.OutputFile(alignment_metrics_tar),
                lib,
            ),
            kwargs={'realign': args['realign']}
        )

    if run_hmmcopy:
        if not run_alignment:
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
                mgd.OutputFile(reads_csvs),
                mgd.OutputFile(segs_csvs),
                mgd.OutputFile(metrics_csvs),
                mgd.OutputFile(params_csvs),
                mgd.OutputFile(igv_csvs),
                mgd.OutputFile(segs_pdf),
                mgd.OutputFile(bias_pdf),
                mgd.OutputFile(heatmap_pdf),
                mgd.OutputFile(heatmap_filt_pdf),
                mgd.OutputFile(metrics_pdf),
                mgd.OutputFile(kernel_density_pdf),
                mgd.OutputFile(hmmcopy_data_tar),
                cellids,
                config['hmmcopy'],
                sampleinfo
            ),
        )

    if run_annotation:
        workflow.subworkflow(
            name='annotation_workflow',
            ctx={'docker_image': config['annotation']['docker']['single_cell_pipeline']},
            func=qc_annotation.create_qc_annotation_workflow,
            args=(
                mgd.InputFile(metrics_csvs),
                mgd.InputFile(reads_csvs),
                mgd.InputFile(alignment_metrics_csv),
                mgd.InputFile(gc_metrics_csv),
                mgd.InputFile(segs_pdf),
                mgd.OutputFile(merged_metrics_csvs),
                mgd.OutputFile(qc_report),
                mgd.OutputFile(corrupt_tree_newick),
                mgd.OutputFile(consensus_tree_newick),
                mgd.OutputFile(phylo_csv),
                mgd.OutputFile(loci_rank_trees),
                mgd.OutputFile(filtered_data),
                mgd.OutputFile(corrupt_tree_pdf),
                mgd.OutputFile(segs_pass),
                mgd.OutputFile(segs_fail),
                mgd.OutputFile(corrupt_heatmap_pdf),
                config['annotation'],
                lib,
            ),
            kwargs={'no_corrupt_tree': args['no_corrupt_tree']}
        )

    return workflow
