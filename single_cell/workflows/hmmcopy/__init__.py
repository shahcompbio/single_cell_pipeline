'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import copy
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers

def create_hmmcopy_workflow(
        bam_file, bai_file, hmmcopy_data, igv_seg_filename, segs_pdf, bias_pdf,
        plot_heatmap_ec_output,
        plot_heatmap_ec_filt_output, plot_metrics_output,
        plot_kernel_density_output, meta_yaml, cell_ids, config, args,
        hmmparams, params_tag, results_dir, alignment_metrics=None):

    sample_info = helpers.get_sample_info(args["input_yaml"])

    chromosomes = config["chromosomes"]

    multipliers = copy.deepcopy(hmmparams["multipliers"])
    multipliers.append(0)

    rows = [int(cellinfo["row"]) for cellinfo in sample_info.values()]
    rows = sorted(set(rows))

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cell_ids,
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('sampleinfo', 'cell_id', axes_origin=[]),
        value=sample_info)

    workflow.setobj(
        obj=mgd.OutputChunks('row'),
        value=rows,
    )


    workflow.transform(
        name='run_hmmcopy',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.hmmcopy.tasks.run_hmmcopy",
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
            mgd.TempSpace('hmmcopy_temp', 'cell_id'),
            mgd.TempInputObj('sampleinfo', 'cell_id'),
        ),
    )

    workflow.transform(
        name='merge_reads',
        ctx={
            'mem': config["memory"]['low'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.hmmcopy.tasks.merge_hdf_files_on_disk",
        args=(
            mgd.TempInputFile('reads.h5', 'cell_id'),
            mgd.TempOutputFile("reads.h5"),
            multipliers,
            'hmmcopy/reads'
        ),
        kwargs={
            'dtypes': {'valid': bool, 'ideal': bool, 'state':float, 'multiplier':float}
        }
    )

    workflow.transform(
        name='merge_segs',
        ctx={
            'mem': config["memory"]['low'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.hmmcopy.tasks.merge_hdf_files_on_disk",
        args=(
            mgd.TempInputFile('segs.h5', 'cell_id'),
            mgd.TempOutputFile("segments.h5"),
            multipliers,
            'hmmcopy/segments'
        ),
        kwargs={
            'dtypes': {'end': int, 'median': float}
        }
    )

    workflow.transform(
        name='merge_metrics',
        ctx={
            'mem': config["memory"]['low'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.hmmcopy.tasks.merge_hdf_files_in_memory",
        args=(
            mgd.TempInputFile('hmm_metrics.h5', 'cell_id'),
            mgd.TempOutputFile("hmmcopy_metrics.h5"),
            multipliers,
            'hmmcopy/metrics'
        ),
        kwargs={
            'dtypes': {'too_even': float,
                       'mad_neutral_state': float,
                       }
        }
    )

    workflow.transform(
        name='merge_params',
        ctx={
            'mem': config["memory"]['low'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.hmmcopy.tasks.merge_hdf_files_in_memory",
        args=(
            mgd.TempInputFile('params.h5', 'cell_id'),
            mgd.TempOutputFile("params.h5"),
            multipliers,
            'hmmcopy/params'
        ),
    )

    annotation_input = 'hmmcopy_metrics.h5'
    if alignment_metrics:
        annotation_input = 'hmmcopy_quality_metrics.h5'
        workflow.transform(
            name="add_quality",
            ctx={
                'mem': config["memory"]['low'],
                'pool_id': config['pools']['standard'],
                'ncpus': 1},
            func="single_cell.workflows.hmmcopy.tasks.add_quality",
            args=(
                mgd.TempInputFile('hmmcopy_metrics.h5'),
                mgd.InputFile(alignment_metrics),
                multipliers,
                mgd.TempOutputFile("hmmcopy_quality_metrics.h5"),
                hmmparams['classifier_training_data'],
            ),
        )

    workflow.transform(
        name='annotate_metrics_with_info_and_clustering',
        ctx={
            'mem': config["memory"]['low'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.hmmcopy.tasks.annotate_metrics",
        args=(
            mgd.TempInputFile('reads.h5'),
            mgd.TempInputFile(annotation_input),
            mgd.TempOutputFile("annotated_metrics.h5"),
            sample_info,
            cell_ids,
            multipliers,
        ),
        kwargs={'chromosomes': config["chromosomes"]}
    )

    workflow.transform(
        name='merge_hmm_copy_plots',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.hmmcopy.tasks.merge_pdf",
        args=(
            [
                mgd.TempInputFile('segments.pdf', 'cell_id'),
                mgd.TempInputFile('bias.pdf', 'cell_id'),
            ],
            [
                mgd.OutputFile("segs", "row", axes_origin=[], template=segs_pdf),
                mgd.OutputFile("bias", "row", axes_origin=[], template=bias_pdf),
            ],
            mgd.TempInputFile("annotated_metrics.h5"),
            None,
            rows
        )
    )

    workflow.transform(
        name='create_igv_seg',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.hmmcopy.tasks.create_igv_seg",
        args=(
            mgd.TempInputFile("segments.h5"),
            mgd.TempInputFile("annotated_metrics.h5"),
            mgd.OutputFile(igv_seg_filename),
            hmmparams
        )
    )

    workflow.transform(
        name='plot_metrics',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.hmmcopy.tasks.plot_metrics",
        args=(
            mgd.TempInputFile("annotated_metrics.h5"),
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
        func="single_cell.workflows.hmmcopy.tasks.plot_kernel_density",
        args=(
            mgd.TempInputFile('annotated_metrics.h5'),
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
        func="single_cell.workflows.hmmcopy.tasks.plot_pcolor",
        args=(
            mgd.TempInputFile('reads.h5'),
            mgd.TempInputFile('annotated_metrics.h5'),
            mgd.OutputFile(plot_heatmap_ec_output),
            mgd.TempSpace("heatmap_ec_temp"),
            multipliers
        ),
        kwargs={
            'plot_title': 'QC pipeline metrics',
            'column_name': 'state',
            'plot_by_col': 'experimental_condition',
            'color_by_col': 'cell_call',
            'chromosomes': chromosomes,
            'max_cn': hmmparams['num_states'],
            'scale_by_cells': False,
            'mappability_threshold': hmmparams["map_cutoff"]
        }
    )

    workflow.transform(
        name='plot_heatmap_ec_filtered',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.hmmcopy.tasks.plot_pcolor",
        args=(
            mgd.TempInputFile('reads.h5'),
            mgd.TempInputFile('annotated_metrics.h5'),
            mgd.OutputFile(plot_heatmap_ec_filt_output),
            mgd.TempSpace("heatmap_ec_filt_temp"),
            multipliers
        ),
        kwargs={
            'plot_title': 'QC pipeline metrics',
            'column_name': 'state',
            'plot_by_col': 'experimental_condition',
            'color_by_col': 'cell_call',
            'chromosomes': chromosomes,
            'max_cn': hmmparams['num_states'],
            'scale_by_cells': False,
            'cell_filters': config["good_cells"],
            'mappability_threshold': hmmparams["map_cutoff"]
        }
    )

    workflow.transform(
        name='merge_all_hdf5_stores',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.hmmcopy.tasks.merge_tables",
        args=(
            mgd.TempInputFile("reads.h5"),
            mgd.TempInputFile("segments.h5"),
            mgd.TempInputFile("annotated_metrics.h5"),
            mgd.TempInputFile("params.h5"),
            mgd.OutputFile(hmmcopy_data),
        )
    )

    metadata = {
        'hmmcopy':{
            'data': hmmcopy_data,
            'reads_table': '/hmmcopy/reads/0',
            'parameters_table': '/hmmcopy/params/0',
            'segments_table': '/hmmcopy/segments/0',
            'metrics_table': '/hmmcopy/metrics/0',
            'hmmcopy_params_tag': params_tag,
            'hmmcopy_params': hmmparams,
            'chromosomes': config['chromosomes'],
            'ref_genome': config['ref_genome'],
            'cell_filters': config["good_cells"]
        }
    }


    workflow.transform(
        name='generate_meta_yaml',
        ctx={
            'mem': config["memory"]['low'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.utils.helpers.write_to_yaml",
        args=(
            mgd.OutputFile(meta_yaml),
            metadata
        )
    )


    return workflow
