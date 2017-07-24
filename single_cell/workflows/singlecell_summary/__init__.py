'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
import tasks


def create_summary_workflow(hmm_segments, hmm_reads, hmm_metrics, metrics_summary, gc_matrix, cn_matrix, config, args, sample_ids):


    results_dir = os.path.join(args['out_dir'], 'results')

    hmmcopy_segments_filename = os.path.join(results_dir, 'segments.csv')
    hmmcopy_reads_filename = os.path.join(results_dir, 'reads.csv')
    hmmcopy_hmm_reads_filt_filename = os.path.join(results_dir, 'filtered_reads.csv')
    hmmcopy_hmm_segs_filt_filename = os.path.join(results_dir, 'filtered_segs.csv')
    all_metrics_heatmap_filename = os.path.join(results_dir, 'all_metrics_summary_hmap.csv')
    gc_metrics_filename = os.path.join(results_dir, 'gc_metrics_summary.csv')
    cn_matrix_filename = os.path.join(results_dir, 'cn_matrix.csv')

    reads_plot_filename = os.path.join(results_dir, 'plots', 'corrected_reads.pdf')
    bias_plot_filename = os.path.join(results_dir, 'plots', 'bias.pdf')
    segs_plot_filename = os.path.join(results_dir, 'plots', 'segments.pdf')

    reads_plot_filename_mad = os.path.join(results_dir, 'plots', 'corrected_reads_mad_0.2.pdf')
    bias_plot_filename_mad = os.path.join(results_dir, 'plots', 'bias_mad_0.2.pdf')
    segs_plot_filename_mad = os.path.join(results_dir, 'plots', 'segments_mad_0.2.pdf')

    order_data_all_output = os.path.join(results_dir, 'plots', 'plot_heatmap_all.csv')

    plot_heatmap_ec_output = os.path.join(results_dir, 'plots', 'plot_heatmap_ec.pdf')
    plot_heatmap_ec_mad_output = os.path.join(results_dir, 'plots', 'plot_heatmap_ec_mad.pdf')
    plot_heatmap_ec_numreads_output = os.path.join(results_dir, 'plots', 'plot_heatmap_ec_numreads.pdf')


    plot_metrics_output = os.path.join(results_dir, 'plots', 'plot_metrics.pdf')
    plot_kernel_density_output = os.path.join(results_dir, 'plots', 'plot_kernel_density.pdf')
    summary_metrics_output = os.path.join(results_dir, 'plots', 'summary_metrics.txt')


    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )
    
    workflow.transform(
        name='merge_tables',
        ctx={'mem': config['med_mem']},
        func=tasks.concatenate_csv,
        args=(
            mgd.InputFile('hmm_segments', 'sample_id', fnames=hmm_segments),
            mgd.OutputFile(hmmcopy_segments_filename),
        ),
    )

    workflow.transform(
        name='merge_reads',
        ctx={'mem': config['high_mem']},
        func=tasks.concatenate_csv,
        args=(
            mgd.InputFile('hmm_reads', 'sample_id', fnames=hmm_reads),
            mgd.OutputFile(hmmcopy_reads_filename),
        ),
    )

    workflow.transform(
        name='merge_hmm_metrics',
        ctx={'mem': config['low_mem']},
        func=tasks.concatenate_csv,
        args=(
            mgd.InputFile('hmm_metrics', 'sample_id', fnames=hmm_metrics),
            mgd.TempOutputFile('hmmcopy_hmm_metrics.csv'),
        ),
    )

    workflow.transform(
        name='merge_summary_metrics',
        ctx={'mem': config['low_mem']},
        func=tasks.concatenate_csv,
        args=(
            mgd.InputFile('metrics_summary', 'sample_id', fnames=metrics_summary),
            mgd.TempOutputFile('metrics_summary.csv'),
        ),
    )

    workflow.transform(
        name='merge_gc_metrics',
        ctx={'mem': config['low_mem']},
        func=tasks.merge_csv,
        args=(
            mgd.InputFile('gc_matrix', 'sample_id', fnames=gc_matrix),
            mgd.OutputFile(gc_metrics_filename),
            'outer',
            'gc'
        ),
    )

    workflow.transform(
        name='merge_cn_metrics',
        ctx={'mem': config['low_mem']},
        func=tasks.merge_csv,
        args=(
            mgd.InputFile('cn_matrix', 'sample_id', fnames=cn_matrix),
            mgd.OutputFile(cn_matrix_filename),
            'outer',
            'chr,start,end,width'
        ),
    )

    workflow.transform(
        name='filter_hmmcopy_results',
        ctx={'mem': config['high_mem']},
        func=tasks.filter_hmm_data,
        args=(
            mgd.TempInputFile('hmmcopy_hmm_metrics.csv'),
            mgd.InputFile(hmmcopy_segments_filename),
            mgd.InputFile(hmmcopy_reads_filename),
            0.2,
            mgd.OutputFile(hmmcopy_hmm_reads_filt_filename),
            mgd.OutputFile(hmmcopy_hmm_segs_filt_filename),
            )
        )

    workflow.transform(
        name='merge_all_metrics',
        ctx={'mem': config['low_mem']},
        func=tasks.merge_tables,
        args=(
            [mgd.TempInputFile('metrics_summary.csv'), mgd.TempInputFile('hmmcopy_hmm_metrics.csv')],
            mgd.TempOutputFile('all_metrics.csv'),
            'merge', ',', 'outer', 'cell_id', 'NA'
            )
        )

    workflow.transform(
        name='plot_heatmap_all',
        ctx={'mem': config['high_mem']},
        func=tasks.plot_heatmap,
        args=(
             mgd.InputFile(hmmcopy_reads_filename),
             mgd.TempInputFile('all_metrics.csv'),
             mgd.OutputFile(order_data_all_output),
             None,
            ),
            kwargs={
                    'plot_title':'QC pipeline metrics',
                    'colname':'integer_copy_number',
                    }

        )

    workflow.transform(
        name='merge_all_metrics_heatmap',
        ctx={'mem': config['high_mem']},
        func=tasks.merge_tables,
        args=(
              [mgd.TempInputFile('all_metrics.csv'), mgd.InputFile(order_data_all_output)],
              mgd.OutputFile(all_metrics_heatmap_filename),
              'merge', ',', 'outer', 'cell_id', 'NA'
            )
        )

    workflow.transform(
        name='plot_hmm_copy',
        ctx={'mem': config['high_mem']},
        func=tasks.plot_hmmcopy,
        args=(
            mgd.InputFile(hmmcopy_reads_filename),
            mgd.InputFile(hmmcopy_segments_filename),
            mgd.InputFile(all_metrics_heatmap_filename),
            mgd.InputFile(config['ref_genome']),
            mgd.OutputFile(reads_plot_filename),
            mgd.OutputFile(segs_plot_filename),
            mgd.OutputFile(bias_plot_filename),
        ),
        kwargs = {
                  'num_states':config['num_states'],
                  'plot_title':'QC pipeline metrics',
                  }

    )

    workflow.transform(
        name='plot_hmm_copy_mad',
        ctx={'mem': config['high_mem']},
        func=tasks.plot_hmmcopy,
        args=(
            mgd.InputFile(hmmcopy_reads_filename),
            mgd.InputFile(hmmcopy_segments_filename),
            mgd.InputFile(all_metrics_heatmap_filename),
            mgd.InputFile(config['ref_genome']),
            mgd.OutputFile(reads_plot_filename_mad),
            mgd.OutputFile(segs_plot_filename_mad),
            mgd.OutputFile(bias_plot_filename_mad),
        ),
        kwargs = {
                  'num_states':config['num_states'],
                  'plot_title':'QC pipeline metrics',
                  'mad_threshold':config['hmmcopy_plot_mad_threshold'],
                  }
    )


    workflow.transform(
        name='plot_metrics',
        ctx={'mem': config['high_mem']},
        func=tasks.plot_metrics,
        args=(
            mgd.InputFile(all_metrics_heatmap_filename),
            mgd.OutputFile(plot_metrics_output),
            'QC pipeline metrics',
            mgd.InputFile(gc_metrics_filename),
            mgd.InputFile(config['gc_windows'])
            )
        )

    workflow.transform(
        name='plot_kernel_density',
        ctx={'mem': config['high_mem']},
        func=tasks.plot_kernel_density,
        args=(
            mgd.InputFile(all_metrics_heatmap_filename),
            mgd.OutputFile(plot_kernel_density_output),
            ',',
            'mad_neutral_state'
            'QC pipeline metrics'
            )
        )

    workflow.transform(
        name='summary_metrics',
        ctx={'mem': config['low_mem']},
        func=tasks.get_summary_metrics,
        args=(
            mgd.InputFile(all_metrics_heatmap_filename),
            mgd.OutputFile(summary_metrics_output),
            )
        )

    workflow.transform(
        name='plot_heatmap_ec',
        ctx={'mem': config['high_mem']},
        func=tasks.plot_heatmap,
        args=(
             mgd.InputFile(hmmcopy_reads_filename),
             mgd.InputFile(all_metrics_heatmap_filename),
             None,
             mgd.OutputFile(plot_heatmap_ec_output),
            ),
            kwargs={
                    'plot_title':'QC pipeline metrics',
                    'colname':'integer_copy_number',
                    'plot_by_col':'experimental_condition',
                    }
        )

    workflow.transform(
        name='plot_heatmap_ec_mad',
        ctx={'mem': config['high_mem']},
        func=tasks.plot_heatmap,
        args=(
             mgd.InputFile(hmmcopy_reads_filename),
             mgd.InputFile(all_metrics_heatmap_filename),
             None,
             mgd.OutputFile(plot_heatmap_ec_mad_output),
            ),
            kwargs={
                    'plot_title':'QC pipeline metrics',
                    'colname':'integer_copy_number',
                    'plot_by_col':'experimental_condition',
                    'mad_threshold':config['heatmap_plot_mad_threshold'],
                    }
        )


    workflow.transform(
        name='plot_heatmap_ec_nreads',
        ctx={'mem': config['high_mem']},
        args=(
             mgd.InputFile(hmmcopy_reads_filename),
             mgd.InputFile(all_metrics_heatmap_filename),
             None,
             mgd.OutputFile(plot_heatmap_ec_numreads_output),
            ),
            kwargs={
                    'plot_title':'QC pipeline metrics',
                    'colname':'integer_copy_number',
                    'plot_by_col':'experimental_condition',
                    'numreads_threshold':config['heatmap_plot_numreads_threshold']
                    }
        )
    return workflow



