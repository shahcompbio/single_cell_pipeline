'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd

import tasks
from single_cell.utils import helpers


def create_copyclone_workflow(bam_file, bai_file, outputdata,
                              segments_plot, bias_plot, plot_metrics_output,
                              plot_kernel_density_output, plot_heatmap_ec_output,
                              plot_heatmap_ec_filt_output,
                              sample_ids, config, args, results_dir):

    sample_info = helpers.get_sample_info(args["input_yaml"])

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )

    workflow.transform(
        name='correct_reads',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func=tasks.correct_reads,
        axes=('sample_id',),
        args=(
            mgd.InputFile('bam_markdups', 'sample_id', fnames=bam_file),
            mgd.InputFile('bai_markdups', 'sample_id', fnames=bai_file),
            mgd.TempOutputFile('corrected_reads.csv', 'sample_id'),
            config,
            mgd.TempSpace('hmmcopy_temp', 'sample_id'),
            mgd.InputInstance("sample_id")
        ),
    )

    workflow.transform(
        name='merge_corrected_reads',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func=tasks.merge_reads,
        args=(
            mgd.TempInputFile('corrected_reads.csv', 'sample_id'),
            mgd.TempOutputFile('corrected_reads.csv'),
        ),
    )

    workflow.transform(
        name="run_copyclone",
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func=tasks.run_copyclone,
        args=(
            mgd.TempInputFile('corrected_reads.csv'),
            mgd.TempOutputFile('reads.h5'),
            mgd.TempOutputFile('segments.h5'),
            mgd.TempOutputFile("metrics.h5"),
            mgd.OutputFile(segments_plot),
            mgd.OutputFile(bias_plot),
            mgd.TempSpace("copyclone_temp"),
            config,
            sample_ids,
            sample_info
        ),
        kwargs={
            "tau": config["copyclone"]["tau"],
            "nu": config["copyclone"]["nu"],
            "eta": config["copyclone"]["eta"],
            "shape": config["copyclone"]["shape"],
            "rate": config["copyclone"]["rate"],
            "ploidy_states": config["copyclone"]["ploidy_states"],
            "num_states": config["copyclone"]["num_states"],
            "num_cores": config["max_cores"]
        }
    )

    workflow.transform(
        name='annotate_metrics_with_info_and_clustering',
        ctx={
            'mem': config["memory"]['low'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.copyclone.tasks.annotate_metrics",
        args=(
            mgd.TempInputFile('reads.h5'),
            mgd.TempInputFile('metrics.h5'),
            mgd.TempOutputFile('metrics_annotated.h5'),
            sample_info,
            sample_ids,
        ),
        kwargs={'chromosomes': config["chromosomes"]}
    )

    workflow.transform(
        name='plot_metrics',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.copyclone.tasks.plot_metrics",
        args=(
            mgd.TempInputFile('metrics_annotated.h5'),
            mgd.OutputFile(plot_metrics_output),
            mgd.TempSpace("plot_metrics_temp"),
            'QC pipeline metrics',
        )
    )

    workflow.transform(
        name='plot_kernel_density',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.copyclone.tasks.plot_kernel_density",
        args=(
            mgd.TempInputFile('metrics_annotated.h5'),
            mgd.OutputFile(plot_kernel_density_output),
            ',',
            'mad_neutral_state',
            'QC pipeline metrics',
        )
    )

    workflow.transform(
        name='plot_heatmap_ec',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.copyclone.tasks.plot_pcolor",
        args=(
            mgd.TempInputFile('reads.h5'),
            mgd.TempInputFile('metrics_annotated.h5'),
            mgd.OutputFile(plot_heatmap_ec_output),
        ),
        kwargs={
            'plot_title': 'QC pipeline metrics',
            'column_name': 'state',
            'plot_by_col': 'experimental_condition',
            'color_by_col': 'cell_call',
            'chromosomes': config['chromosomes'],
            'max_cn': config['copyclone']['num_states'],
            'scale_by_cells': False,
            'mappability_threshold': config['copyclone']["map_cutoff"]
        }
    )

    workflow.transform(
        name='merge_all_hdf5_stores',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.utils.hdfutils.concat_hdf_tables",
        args=(
            [mgd.TempInputFile("reads.h5"),
             mgd.TempInputFile("segments.h5"),
             mgd.TempInputFile("metrics_annotated.h5")],
            mgd.OutputFile(outputdata),
        )
    )

    return workflow
