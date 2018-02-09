'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
import tasks


def create_summary_workflow(sample_info, hmm_segments, hmm_reads, sample_gc_metrics, hmm_metrics,
                            metrics_summary, config, hmmparams, results_dir,
                            args, samples):

    lib = args['library_id']

    gc_metrics = os.path.join(results_dir, '{}_gc_metrics.csv'.format(args['library_id']))

    all_metrics_file = os.path.join(results_dir, '{}_all_metrics_summary.csv'.format(lib))

    plots_dir = os.path.join(results_dir, 'plots')

    plot_heatmap_ec_output = os.path.join(plots_dir, '{}_plot_heatmap_ec.pdf'.format(lib))
    plot_heatmap_ec_mad_output = os.path.join(plots_dir,
                                              '{}_plot_heatmap_ec_mad.pdf'.format(lib))
    plot_heatmap_ec_numreads_output = os.path.join(plots_dir,
                                                   '{}_plot_heatmap_ec_numreads.pdf'.format(lib))

    plot_metrics_output = os.path.join(plots_dir, '{}_plot_metrics.pdf'.format(lib))
    plot_kernel_density_output = os.path.join(plots_dir,
                                              '{}_plot_kernel_density.pdf'.format(lib))
    summary_metrics_output = os.path.join(results_dir, '{}_summary_metrics.txt'.format(lib))

    chromosomes = config["chromosomes"]

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples,
    )


    workflow.transform(
        name='merge_gc_metrics',
        ctx={'mem': config["memory"]['low'], 'pool_id': config['pools']['standard']},
        func=tasks.merge_csv,
        args=(
            mgd.InputFile("gc_metrics.csv", 'sample_id', fnames = sample_gc_metrics),
            mgd.OutputFile(gc_metrics),
            'outer','gc'
        ),
    
    )

    #calculate cell ordering in hierarchical clustering
    workflow.transform(
        name='plot_heatmap_all',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard']},
        func=tasks.plot_heatmap,
        args=(
            mgd.InputFile(hmm_reads),
            None,
            mgd.TempOutputFile('order_data.csv'),
            None,
        ),
        kwargs={
            'plot_title': 'QC pipeline metrics',
            'column_name': 'integer_copy_number',
            'chromosomes': chromosomes,
            'max_cn':hmmparams['num_states'],
        }

    )

    workflow.transform(
        name='merge_all_metrics',
        ctx={'mem': config["memory"]['low'], 'pool_id': config['pools']['standard']},
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
        ctx={'mem': config["memory"]['low'], 'pool_id': config['pools']['standard']},
        func=tasks.plot_metrics,
        args=(
            mgd.InputFile(all_metrics_file),
            mgd.OutputFile(plot_metrics_output),
            'QC pipeline metrics',
            mgd.InputFile(gc_metrics),
            config['gc_windows'],
        )
    )

    workflow.transform(
        name='plot_kernel_density',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard']},
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
        ctx={'mem': config["memory"]['low'], 'pool_id': config['pools']['standard']},
        func=tasks.get_summary_metrics,
        args=(
            mgd.InputFile(all_metrics_file),
            mgd.OutputFile(summary_metrics_output),
        )
    )


    workflow.transform(
        name='plot_heatmap_ec',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard']},
        func=tasks.plot_pcolor,
        args=(
            mgd.InputFile(hmm_reads),
            mgd.InputFile(all_metrics_file),
            None,
            mgd.OutputFile(plot_heatmap_ec_output),
        ),
        kwargs={
            'plot_title': 'QC pipeline metrics',
            'column_name': 'integer_copy_number',
            'plot_by_col': 'experimental_condition',
            'chromosomes': chromosomes,
            'max_cn':hmmparams['num_states'],
        }
    )

    workflow.transform(
        name='plot_heatmap_ec_mad',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard']},
        func=tasks.plot_pcolor,
        args=(
            mgd.InputFile(hmm_reads),
            mgd.InputFile(all_metrics_file),
            None,
            mgd.OutputFile(plot_heatmap_ec_mad_output),
        ),
        kwargs={
            'plot_title': 'QC pipeline metrics',
            'column_name': 'integer_copy_number',
            'plot_by_col': 'experimental_condition',
            'mad_threshold': config['heatmap_plot_mad_threshold'],
            'chromosomes': chromosomes,
            'max_cn':hmmparams['num_states'],
        }
    )

    workflow.transform(
        name='plot_heatmap_ec_nreads',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard']},
        func=tasks.plot_pcolor,
        args=(
            mgd.InputFile(hmm_reads),
            mgd.InputFile(all_metrics_file),
            None,
            mgd.OutputFile(plot_heatmap_ec_numreads_output),
        ),
        kwargs={
            'plot_title': 'QC pipeline metrics',
            'column_name': 'integer_copy_number',
            'plot_by_col': 'experimental_condition',
            'numreads_threshold': config['heatmap_plot_numreads_threshold'],
            'chromosomes': chromosomes,
            'max_cn':hmmparams['num_states'],
        }
    )

    return workflow
