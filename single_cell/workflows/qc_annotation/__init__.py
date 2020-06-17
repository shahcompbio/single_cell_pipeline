'''
Created on Jul 6, 2017

@author: dgrewal
'''

import pypeliner
import pypeliner.managed as mgd


def create_qc_annotation_workflow(
        hmmcopy_metrics, hmmcopy_reads, alignment_metrics, gc_metrics, segs_tar,
        merged_metrics, qc_report, pass_segs, fail_segs,
        plot_heatmap_ec_filt_output, config, alignment_config
):
    ctx = {'docker_image': config['docker']['single_cell_pipeline']}

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.transform(
        name='cell_cycle_classifier',
        func="single_cell.workflows.qc_annotation.tasks.cell_cycle_classifier",
        args=(
            mgd.InputFile(hmmcopy_reads),
            mgd.InputFile(hmmcopy_metrics, extensions=['.yaml']),
            mgd.InputFile(alignment_metrics),
            mgd.TempOutputFile('cell_state_classifier.csv.gz', extensions=['.yaml']),
            mgd.TempSpace('tempdata_cell_cycle')
        ),
        kwargs={'docker_image': config['docker']['cell_cycle_classifier']}
    )

    workflow.transform(
        name="add_quality",
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        func="single_cell.workflows.qc_annotation.tasks.add_quality",
        args=(
            mgd.TempInputFile('cell_state_classifier.csv.gz', extensions=['.yaml']),
            mgd.InputFile(alignment_metrics, extensions=['.yaml']),
            mgd.TempOutputFile("hmmcopy_quality_metrics.csv.gz", extensions=['.yaml']),
            config['classifier_training_data'],
            mgd.TempSpace("hmmcopy_classify_tempdir")
        ),
    )

    workflow.transform(
        name='merge_alignment_hmmcopy_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        func="single_cell.workflows.qc_annotation.tasks.merge_metrics",
        args=(
            mgd.TempInputFile("hmmcopy_quality_metrics.csv.gz", extensions=['.yaml']),
            mgd.InputFile(alignment_metrics, extensions=['.yaml']),
            mgd.TempOutputFile('merged_metrics.csv.gz', extensions=['.yaml'])
        )
    )

    workflow.transform(
        name='add_contamination_status',
        ctx={'mem': config['memory']['med']},
        func="single_cell.workflows.qc_annotation.tasks.add_contamination_status",
        args=(
            mgd.TempInputFile('merged_metrics.csv.gz', extensions=['.yaml']),
            mgd.TempOutputFile('merged_metrics_contamination.csv.gz', extensions=['.yaml']),
            alignment_config['fastq_screen_params']
        ),
        kwargs={
            'reference': config['ref_type'],
        }
    )

    workflow.transform(
        name='generate_qc_report',
        func="single_cell.workflows.qc_annotation.tasks.generate_qc_report",
        args=(
            mgd.TempSpace("QC_report_singlecellpipeline"),
            config['reference_gc'],
            config['fastqscreen_training_data'],
            mgd.TempInputFile('merged_metrics_contamination.csv.gz', extensions=['.yaml']),
            mgd.InputFile(gc_metrics, extensions=['.yaml']),
            mgd.OutputFile(qc_report),
            mgd.TempOutputFile('merged_metrics_contamination_species.csv.gz', extensions=['.yaml']),
        )
    )

    workflow.transform(
        name='filter_segs_plots',
        func="single_cell.workflows.qc_annotation.tasks.filter_plot_tar",
        args=(
            mgd.TempInputFile('merged_metrics_contamination_species.csv.gz', extensions=['.yaml']),
            mgd.InputFile(segs_tar),
            mgd.OutputFile(pass_segs),
            mgd.OutputFile(fail_segs),
            mgd.TempSpace("filter_seg_plots"),
            config['good_cells']
        )
    )

    workflow.transform(
        name='plot_heatmap_ec_filtered',
        func="single_cell.workflows.qc_annotation.tasks.plot_pcolor",
        args=(
            mgd.InputFile(hmmcopy_reads, extensions=['.yaml']),
            mgd.TempInputFile('merged_metrics_contamination_species.csv.gz', extensions=['.yaml']),
            mgd.OutputFile(plot_heatmap_ec_filt_output),
        ),
        kwargs={
            'plot_title': 'QC pipeline metrics',
            'column_name': 'state',
            'plot_by_col': 'experimental_condition',
            'color_by_col': 'cell_call',
            'chromosomes': config['chromosomes'],
            'max_cn': config['num_states'],
            'scale_by_cells': False,
            'cell_filters': config["good_cells"],
            'mappability_threshold': config["map_cutoff"]
        }
    )

    workflow.transform(
        name='finalize_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'num_retry': 1},
        func="single_cell.utils.csvutils.rewrite_csv_file",
        args=(
            mgd.TempInputFile('merged_metrics_contamination_species.csv.gz', extensions=['.yaml']),
            mgd.OutputFile(merged_metrics, extensions=['.yaml']),
        ),
    )

    return workflow
