'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd

import tasks


import copy

from single_cell.utils import helpers


def create_hmmcopy_workflow(
        bam_file, bai_file, hmmcopy_data, igv_seg_filename, segs_pdf, bias_pdf,
        segs_filt_pdf, bias_filt_pdf, plot_heatmap_ec_output,
        plot_heatmap_ec_filt_output, plot_metrics_output,
        plot_kernel_density_output, cell_ids, config, args, hmmparams,
        results_dir):

    sample_info = helpers.get_sample_info(args["input_yaml"])

    chromosomes = config["chromosomes"]

    multipliers = copy.deepcopy(hmmparams["multipliers"])
    multipliers.append(0)

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cell_ids,
    )

    workflow.transform(
        name='run_hmmcopy',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func=tasks.run_hmmcopy,
        axes=('cell_id',),
        args=(
            mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_file),
            mgd.InputFile('bai_markdups', 'cell_id', fnames=bai_file),
            mgd.TempOutputFile('reads.h5', 'cell_id'),
            mgd.TempOutputFile('segs.h5', 'cell_id'),
            mgd.TempOutputFile('params.h5', 'cell_id'),
            mgd.TempOutputFile('hmm_metrics.h5', 'cell_id'),
            mgd.TempOutputFile('segments.pdf', 'cell_id'),
            mgd.TempOutputFile('bias.pdf', 'cell_id'),
            mgd.InputInstance('cell_id'),
            config,
            hmmparams,
            multipliers,
            mgd.TempSpace('hmmcopy_temp', 'cell_id')
        ),
    )

    workflow.transform(
        name='merge_files',
        ctx={
            'mem': config["memory"]['low'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func=tasks.merge_files,
        args=(
            mgd.TempInputFile('reads.h5', 'cell_id'),
            mgd.TempInputFile('segs.h5', 'cell_id'),
            mgd.TempInputFile('hmm_metrics.h5', 'cell_id'),
            mgd.TempInputFile('params.h5', 'cell_id'),
            mgd.TempOutputFile("segments.h5"),
            mgd.TempOutputFile("reads.h5"),
            mgd.TempOutputFile("hmmcopy_metrics.h5"),
            mgd.TempOutputFile("params.h5"),
            mgd.TempSpace("hmmcopy_merge_files_temp"),
            mgd.OutputFile(igv_seg_filename),
            sample_info,
            0.2,
            hmmparams,
            multipliers
        ),
    )

    workflow.transform(
        name='merge_hmm_copy_plots',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func=tasks.merge_pdf,
        args=(
            [
                mgd.TempInputFile('segments.pdf', 'cell_id'),
                mgd.TempInputFile('bias.pdf', 'cell_id'),
            ],
            [
                mgd.OutputFile(segs_pdf),
                mgd.OutputFile(bias_pdf),
            ],
            mgd.TempInputFile("hmmcopy_metrics.h5"),
            None,
            None,
            None,
        )
    )

    workflow.transform(
        name='merge_hmm_copy_plots_filtered',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func=tasks.merge_pdf,
        args=(
            [
                mgd.TempInputFile('segments.pdf', 'cell_id'),
                mgd.TempInputFile('bias.pdf', 'cell_id'),
            ],
            [
                mgd.OutputFile(segs_filt_pdf),
                mgd.OutputFile(bias_filt_pdf),
            ],
            mgd.TempInputFile("hmmcopy_metrics.h5"),
            config['plot_mad_threshold'],
            config['plot_numreads_threshold'],
            config['plot_median_hmmcopy_reads_per_bin_threshold']
        )
    )

    workflow.transform(
        name='add_clustering_order',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func=tasks.add_clustering_order,
        args=(
            mgd.TempInputFile('reads.h5'),
            mgd.TempInputFile("hmmcopy_metrics.h5"),
            mgd.TempOutputFile('hmmcopy_metrics_with_cluster_order.h5'),
            multipliers,
            cell_ids,
        ),
    )

    workflow.transform(
        name='plot_metrics',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func=tasks.plot_metrics,
        args=(
            mgd.TempInputFile("hmmcopy_metrics_with_cluster_order.h5"),
            mgd.OutputFile(plot_metrics_output),
            mgd.TempSpace("plot_metrics_temp"),
            'QC pipeline metrics',
            multipliers
        )
    )

    workflow.transform(
        name='plot_kernel_density',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func=tasks.plot_kernel_density,
        args=(
            mgd.TempInputFile('hmmcopy_metrics_with_cluster_order.h5'),
            mgd.OutputFile(plot_kernel_density_output),
            mgd.TempSpace("hmmcopy_kde_plot_temp"),
            ',',
            'mad_neutral_state',
            'QC pipeline metrics',
            multipliers
        )
    )

    workflow.transform(
        name='plot_heatmap_ec',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func=tasks.plot_pcolor,
        args=(
            mgd.TempInputFile('reads.h5'),
            mgd.TempInputFile('hmmcopy_metrics_with_cluster_order.h5'),
            mgd.OutputFile(plot_heatmap_ec_output),
            mgd.TempSpace("heatmap_ec_temp"),
            multipliers
        ),
        kwargs={
            'plot_title': 'QC pipeline metrics',
            'column_name': 'state',
            'plot_by_col': 'experimental_condition',
            'chromosomes': chromosomes,
            'max_cn': hmmparams['num_states'],
        }
    )

    workflow.transform(
        name='plot_heatmap_ec_filtered',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func=tasks.plot_pcolor,
        args=(
            mgd.TempInputFile('reads.h5'),
            mgd.TempInputFile('hmmcopy_metrics_with_cluster_order.h5'),
            mgd.OutputFile(plot_heatmap_ec_filt_output),
            mgd.TempSpace("heatmap_ec_filt_temp"),
            multipliers
        ),
        kwargs={
            'plot_title': 'QC pipeline metrics',
            'column_name': 'state',
            'plot_by_col': 'experimental_condition',
            'numreads_threshold': config['plot_numreads_threshold'],
            'median_hmmcopy_reads_per_bin_threshold': config['plot_median_hmmcopy_reads_per_bin_threshold'],
            'mad_threshold': config['plot_mad_threshold'],
            'chromosomes': chromosomes,
            'max_cn': hmmparams['num_states'],
        }
    )

    workflow.transform(
        name='merge_all_hdf5_stores',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func=tasks.merge_tables,
        args=(
            mgd.TempInputFile("reads.h5"),
            mgd.TempInputFile("segments.h5"),
            mgd.TempInputFile("hmmcopy_metrics_with_cluster_order.h5"),
            mgd.TempInputFile("params.h5"),
            mgd.OutputFile(hmmcopy_data),
            multipliers,
            cell_ids
        )
    )

    return workflow
