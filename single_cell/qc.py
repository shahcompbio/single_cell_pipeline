'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
from workflows import align, hmmcopy, qc_annotation
from single_cell.utils import helpers
import copy


def qc_workflow(args):

    run_alignment = args['alignment']
    run_hmmcopy = args['hmmcopy']
    run_annotation = args['annotation']

    if not any((run_alignment, run_hmmcopy, run_annotation)):
        run_alignment = True
        run_hmmcopy = True
        run_annotation = True

    config = helpers.load_config(args)

    align_config = config['alignment']

    sampleinfo = helpers.get_sample_info(args['input_yaml'])
    cellids = helpers.get_samples(args['input_yaml'])
    bam_files, _ = helpers.get_bams(args['input_yaml'])

    lib = args["library_id"]

    ctx = {'docker_image': align_config['docker']['single_cell_pipeline']}
    workflow = pypeliner.workflow.Workflow(ctx=ctx)


    outdir = os.path.join(args["out_dir"], "results", "QC")

    alignment_metrics_csv = os.path.join(outdir, '{}_alignment_metrics.csv.gz'.format(lib))
    gc_metrics_csv = os.path.join(outdir, '{}_gc_metrics.csv.gz'.format(lib))
    plots_dir = os.path.join(outdir,  'plots')
    plot_metrics_output = os.path.join(plots_dir, '{}_plot_metrics.pdf'.format(lib))

    baseimage = align_config['docker']['single_cell_pipeline']

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
            ctx={'docker_image': baseimage},
            func=align.create_alignment_workflow,
            args=(
                mgd.InputFile('fastq_1', 'cell_id', 'lane', fnames=fastq1_files, axes_origin=[]),
                mgd.InputFile('fastq_2', 'cell_id', 'lane', fnames=fastq2_files, axes_origin=[]),
                mgd.OutputFile('bam_markdups', 'cell_id', fnames=bam_files, axes_origin=[], extensions=['.bai']),
                mgd.TempOutputFile('biobloom_count_metrics', 'cell_id', axes_origin=[]),
                mgd.OutputFile(alignment_metrics_csv),
                mgd.OutputFile(gc_metrics_csv),
                mgd.OutputFile(plot_metrics_output),
                align_config['ref_genome'],
                align_config,
                args,
                triminfo,
                centerinfo,
                sampleinfo,
                cellids,
            ),
        )

    if run_hmmcopy:

        if not args['alignment']:
            workflow.setobj(
                obj=mgd.OutputChunks('cell_id'),
                value=list(bam_files.keys()),
            )

        hmmcopy_config = config['hmmcopy']

        results_dir = os.path.join(args['out_dir'], 'results', 'hmmcopy_autoploidy')

        reads_csvs = os.path.join(results_dir, '{0}_reads.csv.gz'.format(lib))
        segs_csvs = os.path.join(results_dir, '{0}_segments.csv.gz'.format(lib))
        params_csvs = os.path.join(results_dir, '{0}_params.csv.gz'.format(lib))
        metrics_csvs = os.path.join(results_dir, '{0}_metrics.csv.gz'.format(lib))
        igv_csvs = os.path.join(results_dir, '{0}_igv_segments.seg'.format(lib))

        plots_dir = os.path.join(results_dir, "plots")
        segs_pdf = os.path.join(
            plots_dir, "segments", '{}_segs.tar.gz'.format(lib))
        bias_pdf = os.path.join(plots_dir, "bias", '{}_bias.tar.gz'.format(lib))

        heatmap_filt_pdf = os.path.join(
            plots_dir, '{}_heatmap_by_ec_filtered.pdf'.format(lib))
        heatmap_pdf = os.path.join(
            plots_dir, '{}_heatmap_by_ec.pdf'.format(lib))
        metrics_pdf = os.path.join(
            plots_dir, '{}_metrics.pdf'.format(lib))
        kernel_density_pdf = os.path.join(
            plots_dir, '{}_kernel_density.pdf'.format(lib))

        workflow.subworkflow(
            name='hmmcopy_workflow',
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
                cellids,
                hmmcopy_config,
                sampleinfo
            ),
        )

    if run_annotation:
        results_dir = os.path.join(args['out_dir'], 'results', 'QC')
        metrics_csvs = os.path.join(results_dir, 'hmmcopy_autoploidy', '{0}_metrics.csv.gz'.format(lib))
        merged_metrics_csvs = os.path.join(results_dir, '{0}_metrics.csv.gz'.format(lib))
        qc_report = os.path.join(results_dir, '{0}_QC_report.html'.format(lib))

        workflow.subworkflow(
            name='annotation_workflow',
            func=qc_annotation.create_qc_annotation_workflow,
            args=(
                mgd.InputFile(metrics_csvs),
                mgd.OutputFile(alignment_metrics_csv),
                mgd.OutputFile(gc_metrics_csv),
                mgd.OutputFile(merged_metrics_csvs),
                mgd.OutputFile(qc_report),
                config['annotation'],
            ),
        )

    return workflow
