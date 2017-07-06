'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
import tasks


def create_summary_workflow(hmm_segments, hmm_reads, hmm_metrics, metrics_summary, gc_matrix, cn_matrix, config, args):


    scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
    plot_hmmcopy_script = os.path.join(scripts_directory, 'plot_hmmcopy.py')
    plot_heatmap_script = os.path.join(scripts_directory, 'plot_heatmap.py')
    merge_tables_script = os.path.join(scripts_directory, 'merge.py')
    filter_hmmcopy_script = os.path.join(scripts_directory, 'filter_data.py')
    plot_metrics_script = os.path.join(scripts_directory, 'plot_metrics.py')
    plot_kernel_density_script = os.path.join(scripts_directory, 'plot_kernel_density.py')
    summary_metrics_script = os.path.join(scripts_directory, 'summary_metrics.py')



    metrics_directory = os.path.join(args['out_dir'], 'metrics')
    hmmcopy_directory = os.path.join(args['out_dir'], 'hmmcopy')

    hmmcopy_segments_filename = os.path.join(hmmcopy_directory, 'segments.csv')
    hmmcopy_reads_filename = os.path.join(hmmcopy_directory, 'reads.csv')
    hmmcopy_hmm_metrics_filename = os.path.join(hmmcopy_directory, 'hmm_metrics.csv')
    hmmcopy_hmm_reads_filt_filename = os.path.join(hmmcopy_directory, 'filtered_reads.csv')
    hmmcopy_hmm_segs_filt_filename = os.path.join(hmmcopy_directory, 'filtered_segs.csv')

    metrics_summary_filename = os.path.join(metrics_directory, 'summary', 'summary.csv')

    plots_directory = os.path.join(args['out_dir'], 'plots')
    reads_plot_filename = os.path.join(plots_directory, 'corrected_reads.pdf')
    bias_plot_filename = os.path.join(plots_directory, 'bias.pdf')
    segs_plot_filename = os.path.join(plots_directory, 'segments.pdf')

    reads_plot_filename_mad = os.path.join(plots_directory, 'corrected_reads_mad_0.2.pdf')
    bias_plot_filename_mad = os.path.join(plots_directory, 'bias_mad_0.2.pdf')
    segs_plot_filename_mad = os.path.join(plots_directory, 'segments_mad_0.2.pdf')

    plot_heatmap_all_output = os.path.join(plots_directory, 'plot_heatmap_all.pdf')
    order_data_all_output = os.path.join(plots_directory, 'plot_heatmap_all.csv')

    gc_metrics_filename = os.path.join(metrics_directory, 'summary', 'gc_metrics_summary.csv')
    all_metrics_filename = os.path.join(metrics_directory, 'summary', 'all_metrics_summary.csv')
    all_metrics_heatmap_filename = os.path.join(metrics_directory, 'summary', 'all_metrics_summary_hmap.csv')


    plot_heatmap_ec_output = os.path.join(plots_directory, 'plot_heatmap_ec.pdf')
    plot_heatmap_ec_mad_output = os.path.join(plots_directory, 'plot_heatmap_ec_mad.pdf')
    plot_heatmap_ec_numreads_output = os.path.join(plots_directory, 'plot_heatmap_ec_numreads.pdf')

    plot_heatmap_st_output = os.path.join(plots_directory, 'plot_heatmap_st.pdf')
    plot_heatmap_st_mad_output = os.path.join(plots_directory, 'plot_heatmap_st_mad.pdf')
    plot_heatmap_st_numreads_output = os.path.join(plots_directory, 'plot_heatmap_st_numreads.pdf')


    plot_metrics_output = os.path.join(plots_directory, 'plot_metrics.pdf')
    plot_kernel_density_output = os.path.join(plots_directory, 'plot_kernel_density.pdf')
    summary_metrics_output = os.path.join(plots_directory, 'summary_metrics.txt')


    cn_metrics_filename = os.path.join(metrics_directory, 'summary', 'cn_metrics_summary.csv')

    workflow = pypeliner.workflow.Workflow()

    # raise Exception(hmmcopy_segments.values())

    workflow.transform(
        name='merge_tables',
        func=tasks.concatenate_csv,
        args=(
            hmm_segments,
            mgd.OutputFile(hmmcopy_segments_filename),
        ),
    )

    workflow.transform(
        name='merge_reads',
        func=tasks.concatenate_csv,
        args=(
            hmm_reads,
            mgd.OutputFile(hmmcopy_reads_filename),
        ),
    )

    workflow.transform(
        name='merge_hmm_metrics',
        func=tasks.concatenate_csv,
        args=(
            hmm_metrics,
            mgd.OutputFile(hmmcopy_hmm_metrics_filename),
        ),
    )

    workflow.transform(
        name='merge_summary_metrics',
        func=tasks.concatenate_csv,
        args=(
            metrics_summary,
            mgd.OutputFile(metrics_summary_filename),
        ),
    )

    workflow.transform(
        name='merge_gc_metrics',
        func=tasks.merge_csv,
        args=(
            gc_matrix,
            mgd.OutputFile(gc_metrics_filename),
            'outer',
            'gc'
        ),
    )

    workflow.transform(
        name='merge_cn_metrics',
        func=tasks.merge_csv,
        args=(
            cn_matrix,
            mgd.OutputFile(cn_metrics_filename),
            'outer',
            'chr,start,end,width'
        ),
    )

    workflow.commandline(
        name='filter_hmmcopy_results',
        args=(
            config['python'],
            filter_hmmcopy_script,
            '--corrected_reads', mgd.InputFile(hmmcopy_reads_filename),
            '--segments', mgd.InputFile(hmmcopy_segments_filename),
            '--quality_metrics', mgd.InputFile(hmmcopy_hmm_metrics_filename),
            '--reads_output', mgd.OutputFile(hmmcopy_hmm_reads_filt_filename),
            '--segs_output', mgd.OutputFile(hmmcopy_hmm_segs_filt_filename),
            '--mad_threshold', '0.2'
            )
        )



    workflow.commandline(
        name='merge_all_metrics',
        args=(
            config['python'],
            merge_tables_script,
            '--merge_type', 'outer',
            '--nan_value', 'NA',
            '--input', mgd.InputFile(metrics_summary_filename), mgd.InputFile(hmmcopy_hmm_metrics_filename),
            '--key_cols', 'cell_id',
            '--separator', 'comma',
            '--type', 'merge',
            '--output', mgd.OutputFile(all_metrics_filename),
            )
        )

    workflow.commandline(
        name='plot_heatmap_all',
        args=(
            config['python'],
            plot_heatmap_script,
            '--output', mgd.OutputFile(plot_heatmap_all_output),
            '--metrics', mgd.InputFile(all_metrics_filename),
            '--separator', 'comma',
            '--plot_title', 'QC pipeline metrics',
            '--input', mgd.InputFile(hmmcopy_reads_filename),
            '--order_data', mgd.OutputFile(order_data_all_output),
            '--column_name', 'integer_copy_number'
            )
        )

    workflow.commandline(
        name='merge_all_metrics_heatmap',
        args=(
            config['python'],
            merge_tables_script,
            '--merge_type', 'outer',
            '--nan_value', 'NA',
            '--input', mgd.InputFile(all_metrics_filename), mgd.InputFile(order_data_all_output),
            '--key_cols', 'cell_id',
            '--separator', 'comma',
            '--type', 'merge',
            '--output', mgd.OutputFile(all_metrics_heatmap_filename),
            )
        )

    workflow.commandline(
        name='plot_hmm_copy',
        args=(
            config['python'],
            plot_hmmcopy_script,
            '--corrected_reads', mgd.InputFile(hmmcopy_reads_filename),
            '--segments', mgd.InputFile(hmmcopy_segments_filename),
            '--quality_metrics', mgd.InputFile(all_metrics_heatmap_filename),
            '--ref_genome', mgd.InputFile(config['ref_genome']),
            '--num_states', config['num_states'],
            '--reads_output', mgd.OutputFile(reads_plot_filename),
            '--bias_output', mgd.OutputFile(bias_plot_filename),
            '--segs_output', mgd.OutputFile(segs_plot_filename),
            '--plot_title', 'QC pipeline metrics',
        ),
    )

    workflow.commandline(
        name='plot_hmm_copy_mad',
        args=(
            config['python'],
            plot_hmmcopy_script,
            '--corrected_reads', mgd.InputFile(hmmcopy_reads_filename),
            '--segments', mgd.InputFile(hmmcopy_segments_filename),
            '--quality_metrics', mgd.InputFile(all_metrics_heatmap_filename),
            '--ref_genome', mgd.InputFile(config['ref_genome']),
            '--num_states', config['num_states'],
            '--reads_output', mgd.OutputFile(reads_plot_filename_mad),
            '--bias_output', mgd.OutputFile(bias_plot_filename_mad),
            '--segs_output', mgd.OutputFile(segs_plot_filename_mad),
            '--plot_title', 'QC pipeline metrics',
            '--mad_threshold', config['hmmcopy_plot_mad_threshold'],
        ),
    )



    workflow.commandline(
        name='plot_metrics',
        args=(
            config['python'],
            plot_metrics_script,
            mgd.InputFile(all_metrics_filename),
            mgd.OutputFile(plot_metrics_output),
            '--plot_title', 'QC pipeline metrics',
            '--gcbias_matrix', mgd.InputFile(gc_metrics_filename),
            '--gc_content_data', mgd.InputFile(config['gc_windows'])
            )
        )

    workflow.commandline(
        name='plot_kernel_density',
        args=(
            config['python'],
            plot_kernel_density_script,
            '--separator', 'comma',
            '--input', mgd.InputFile(all_metrics_filename),
            '--plot_title', 'QC pipeline metrics',
            '--output', mgd.OutputFile(plot_kernel_density_output),
            '--column_name', 'mad_neutral_state'
            )
        )

    workflow.commandline(
        name='summary_metrics',
        args=(
            config['python'],
            summary_metrics_script,
            '--input', mgd.InputFile(all_metrics_filename),
            '--summary_metrics', mgd.OutputFile(summary_metrics_output),
            )
        )

    workflow.commandline(
        name='plot_heatmap_ec',
        args=(
            config['python'],
            plot_heatmap_script,
            '--output', mgd.OutputFile(plot_heatmap_ec_output),
            '--metrics', mgd.InputFile(all_metrics_filename),
            '--separator', 'comma',
            '--plot_title', 'QC pipeline metrics',
            '--input', mgd.InputFile(hmmcopy_reads_filename),
            '--column_name', 'integer_copy_number',
            '--plot_by_col','experimental_condition',
            )
        )

    workflow.commandline(
        name='plot_heatmap_ec_mad',
        args=(
            config['python'],
            plot_heatmap_script,
            '--output', mgd.OutputFile(plot_heatmap_ec_mad_output),
            '--metrics', mgd.InputFile(all_metrics_filename),
            '--separator', 'comma',
            '--plot_title', 'QC pipeline metrics',
            '--input', mgd.InputFile(hmmcopy_reads_filename),
            '--column_name', 'integer_copy_number',
            '--plot_by_col','experimental_condition',
            '--mad_threshold', config['heatmap_plot_mad_threshold']
            )
        )


    workflow.commandline(
        name='plot_heatmap_ec_nreads',
        args=(
            config['python'],
            plot_heatmap_script,
            '--output', mgd.OutputFile(plot_heatmap_ec_numreads_output),
            '--metrics', mgd.InputFile(all_metrics_filename),
            '--separator', 'comma',
            '--plot_title', 'QC pipeline metrics',
            '--input', mgd.InputFile(hmmcopy_reads_filename),
            '--column_name', 'integer_copy_number',
            '--plot_by_col','experimental_condition',
            '--numreads_threshold', config['heatmap_plot_numreads_threshold']
            )
        )



    workflow.commandline(
        name='plot_heatmap_st',
        args=(
            config['python'],
            plot_heatmap_script,
            '--output', mgd.OutputFile(plot_heatmap_st_output),
            '--metrics', mgd.InputFile(all_metrics_filename),
            '--separator', 'comma',
            '--plot_title', 'QC pipeline metrics',
            '--input', mgd.InputFile(hmmcopy_reads_filename),
            '--column_name', 'integer_copy_number',
            '--plot_by_col','sample_type',
            )
        )

    workflow.commandline(
        name='plot_heatmap_st_mad',
        args=(
            config['python'],
            plot_heatmap_script,
            '--output', mgd.OutputFile(plot_heatmap_st_mad_output),
            '--metrics', mgd.InputFile(all_metrics_filename),
            '--separator', 'comma',
            '--plot_title', 'QC pipeline metrics',
            '--input', mgd.InputFile(hmmcopy_reads_filename),
            '--column_name', 'integer_copy_number',
            '--plot_by_col','sample_type',
            '--mad_threshold', config['heatmap_plot_mad_threshold']
            )
        )


    workflow.commandline(
        name='plot_heatmap_st_nreads',
        args=(
            config['python'],
            plot_heatmap_script,
            '--output', mgd.OutputFile(plot_heatmap_st_numreads_output),
            '--metrics', mgd.InputFile(all_metrics_filename),
            '--separator', 'comma',
            '--plot_title', 'QC pipeline metrics',
            '--input', mgd.InputFile(hmmcopy_reads_filename),
            '--column_name', 'integer_copy_number',
            '--plot_by_col','sample_type',
            '--numreads_threshold', config['heatmap_plot_numreads_threshold']
            )
        )


    return workflow



