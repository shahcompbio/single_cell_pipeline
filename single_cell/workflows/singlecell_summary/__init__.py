'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
import tasks


def create_summary_workflow(sample_info, hmm_segments, hmm_reads, hmm_metrics,
                            metrics_summary, gc_matrix, config,
                            out_dir, sample_ids):

    results_dir = os.path.join(out_dir, 'results')


    all_metrics_file = os.path.join(results_dir, 'all_metrics_summary.csv')

    plots_dir = os.path.join(results_dir, 'plots')

    plot_heatmap_ec_output = os.path.join(plots_dir, 'plot_heatmap_ec.pdf')
    plot_heatmap_ec_mad_output = os.path.join(plots_dir,
                                              'plot_heatmap_ec_mad.pdf')
    plot_heatmap_ec_numreads_output = os.path.join(plots_dir,
                                                   'plot_heatmap_ec_numreads.pdf')

    plot_metrics_output = os.path.join(plots_dir, 'plot_metrics.pdf')
    plot_kernel_density_output = os.path.join(plots_dir,
                                              'plot_kernel_density.pdf')
    summary_metrics_output = os.path.join(plots_dir, 'summary_metrics.txt')

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )

    #calculate cell ordering in hierarchical clustering
    workflow.transform(
        name='plot_heatmap_all',
        ctx={'mem': config['high_mem']},
        func=tasks.plot_heatmap,
        args=(
            mgd.InputFile(hmm_reads),
            None,
            mgd.TempOutputFile('order_data.csv'),
            None,
        ),
        kwargs={
            'plot_title': 'QC pipeline metrics',
            'colname': 'integer_copy_number',
        }

    )

    workflow.transform(
        name='merge_all_metrics',
        ctx={'mem': config['low_mem']},
        func=tasks.merge_tables,
        args=(
            [mgd.InputFile(metrics_summary),
             mgd.InputFile(hmm_metrics),
             mgd.InputFile(sample_info),
             mgd.TempInputFile('order_data.csv')],
            mgd.OutputFile(all_metrics_file),
            'merge', ',', 'outer', 'cell_id', 'NA'
        )
    )

    workflow.transform(
        name='plot_metrics',
        ctx={'mem': config['high_mem']},
        func=tasks.plot_metrics,
        args=(
            mgd.InputFile(all_metrics_file),
            mgd.OutputFile(plot_metrics_output),
            'QC pipeline metrics',
            mgd.InputFile(gc_matrix),
            mgd.InputFile(config['gc_windows'])
        )
    )

    workflow.transform(
        name='plot_kernel_density',
        ctx={'mem': config['high_mem']},
        func=tasks.plot_kernel_density,
        args=(
            mgd.InputFile(all_metrics_file),
            mgd.OutputFile(plot_kernel_density_output),
            ',',
            'mad_neutral_state',
            'QC pipeline metrics'
        )
    )

    workflow.transform(
        name='summary_metrics',
        ctx={'mem': config['low_mem']},
        func=tasks.get_summary_metrics,
        args=(
            mgd.InputFile(all_metrics_file),
            mgd.OutputFile(summary_metrics_output),
        )
    )

    workflow.transform(
        name='plot_heatmap_ec',
        ctx={'mem': config['high_mem']},
        func=tasks.plot_heatmap,
        args=(
            mgd.InputFile(hmm_reads),
            mgd.InputFile(all_metrics_file),
            None,
            mgd.OutputFile(plot_heatmap_ec_output),
        ),
        kwargs={
            'plot_title': 'QC pipeline metrics',
            'colname': 'integer_copy_number',
            'plot_by_col': 'experimental_condition',
        }
    )

    workflow.transform(
        name='plot_heatmap_ec_mad',
        ctx={'mem': config['high_mem']},
        func=tasks.plot_heatmap,
        args=(
            mgd.InputFile(hmm_reads),
            mgd.InputFile(all_metrics_file),
            None,
            mgd.OutputFile(plot_heatmap_ec_mad_output),
        ),
        kwargs={
            'plot_title': 'QC pipeline metrics',
            'colname': 'integer_copy_number',
            'plot_by_col': 'experimental_condition',
            'mad_threshold': config['heatmap_plot_mad_threshold'],
        }
    )

    workflow.transform(
        name='plot_heatmap_ec_nreads',
        ctx={'mem': config['high_mem']},
        func=tasks.plot_heatmap,
        args=(
            mgd.InputFile(hmm_reads),
            mgd.InputFile(all_metrics_file),
            None,
            mgd.OutputFile(plot_heatmap_ec_numreads_output),
        ),
        kwargs={
            'plot_title': 'QC pipeline metrics',
            'colname': 'integer_copy_number',
            'plot_by_col': 'experimental_condition',
            'numreads_threshold': config['heatmap_plot_numreads_threshold']
        }
    )
    return workflow
