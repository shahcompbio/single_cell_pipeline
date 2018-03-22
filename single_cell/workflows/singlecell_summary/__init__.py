'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
import tasks
from single_cell.utils import helpers


def create_summary_workflow(alignment_metrics, gc_metrics, hmm_segments, hmm_reads,
                            hmm_metrics, hmm_params, reads_pdf_output, segs_pdf_output, bias_pdf_output,
                            params_pdf_output, config, hmmparams, results_dir,
                            args, samples):

    lib = args['library_id']

    all_metrics_file = os.path.join(results_dir, '{}_all_metrics_summary.csv'.format(lib))

    plots_dir = os.path.join(results_dir, 'plots')

    plot_heatmap_ec_output = os.path.join(plots_dir, '{}_plot_heatmap_ec.pdf'.format(lib))
    plot_heatmap_ec_filt_output = os.path.join(plots_dir,
                                                   '{}_plot_heatmap_ec_filtered.pdf'.format(lib))

    plot_metrics_output = os.path.join(plots_dir, '{}_plot_metrics.pdf'.format(lib))
    plot_kernel_density_output = os.path.join(plots_dir,
                                              '{}_plot_kernel_density.pdf'.format(lib))
    summary_metrics_output = os.path.join(results_dir, '{}_summary_metrics.txt'.format(lib))

    chromosomes = config["chromosomes"]

    sample_info  = helpers.get_sample_info(args["input_yaml"])

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples,
    )


    #calculate cell ordering in hierarchical clustering
    workflow.transform(
        name='plot_heatmap_all',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
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
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.merge_tables,
        args=(
            [mgd.InputFile(alignment_metrics),
             mgd.InputFile(hmm_metrics),
             mgd.TempInputFile('order_data.csv')],
            mgd.TempOutputFile("all_metrics.csv"),
            ',', 'outer', 'cell_id', 'NA'
        )
    )

    workflow.transform(
        name='annotate_metrics',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.annotate_metrics,
        args=(
            mgd.TempInputFile("all_metrics.csv"),
            sample_info,
            mgd.TempOutputFile("all_metrics_annotated.csv"),
        )
    )

    workflow.transform(
        name='classify_cells',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.classify,
        args=(
            mgd.TempInputFile("all_metrics_annotated.csv"),
            mgd.OutputFile(all_metrics_file),
        )
    )
 
    workflow.transform(
        name='plot_metrics',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
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
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
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
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.get_summary_metrics,
        args=(
            mgd.InputFile(all_metrics_file),
            mgd.OutputFile(summary_metrics_output),
        )
    )
 
 
    workflow.transform(
        name='plot_heatmap_ec',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
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
        name='plot_heatmap_ec_filtered',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.plot_pcolor,
        args=(
            mgd.InputFile(hmm_reads),
            mgd.InputFile(all_metrics_file),
            None,
            mgd.OutputFile(plot_heatmap_ec_filt_output),
        ),
        kwargs={
            'plot_title': 'QC pipeline metrics',
            'column_name': 'integer_copy_number',
            'plot_by_col': 'experimental_condition',
            'numreads_threshold': config['plot_numreads_threshold'],
            'quality_threshold': config['plot_quality_threshold'],
            'mad_threshold': config['plot_mad_threshold'],
            'chromosomes': chromosomes,
            'max_cn':hmmparams['num_states'],
        }
    )


    workflow.transform(
        name='plot_hmm_copy',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.plot_hmmcopy,
        axes=('sample_id',),
        args=(
            mgd.InputFile(hmm_reads),
            mgd.InputFile(hmm_segments),
            mgd.InputFile(hmm_params),
            mgd.InputFile(all_metrics_file),
            None,
            config['ref_genome'],
            mgd.TempOutputFile('reads.pdf', 'sample_id'),
            mgd.TempOutputFile('segs.pdf', 'sample_id'),
            mgd.TempOutputFile('bias.pdf', 'sample_id'),
            mgd.TempOutputFile('params.pdf', 'sample_id'),
            mgd.InputInstance('sample_id'),
        ),
        kwargs={
            'num_states': hmmparams['num_states'],
            'plot_title': 'QC pipeline metrics',
            'annotation_cols': ['cell_call','experimental_condition', 'sample_type',
                                'coverage_depth', 'mad_neutral_state', 'MSRSI_non_integerness'],
        }
    )

    workflow.transform(
        name='merge_hmm_copy',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.merge_pdf,
        args=(
              [mgd.TempInputFile('reads.pdf', 'sample_id'),
              mgd.TempInputFile('segs.pdf', 'sample_id'),
              mgd.TempInputFile('bias.pdf', 'sample_id'),
              mgd.TempInputFile('params.pdf', 'sample_id')],
              [mgd.OutputFile(reads_pdf_output),
              mgd.OutputFile(segs_pdf_output),
              mgd.OutputFile(bias_pdf_output),
              mgd.OutputFile(params_pdf_output)],
              mgd.InputFile(all_metrics_file),
              None,
              None,
              None,
            )
    )

    reads_mad_pdf_output = os.path.join(results_dir, 'plots', '{}_reads_filtered.pdf'.format(lib))
    segs_mad_pdf_output = os.path.join(results_dir, 'plots', '{}_segs_filtered.pdf'.format(lib))
    bias_mad_pdf_output = os.path.join(results_dir, 'plots', '{}_bias_filtered.pdf'.format(lib))
    params_mad_pdf_output = os.path.join(results_dir, 'plots', '{}_params_filtered.pdf'.format(lib))
    workflow.transform(
        name='merge_hmm_copy_filtered',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.merge_pdf,
        args=(
              [mgd.TempInputFile('reads.pdf', 'sample_id'),
              mgd.TempInputFile('segs.pdf', 'sample_id'),
              mgd.TempInputFile('bias.pdf', 'sample_id'),
              mgd.TempInputFile('params.pdf', 'sample_id')],
              [mgd.OutputFile(reads_mad_pdf_output),
              mgd.OutputFile(segs_mad_pdf_output),
              mgd.OutputFile(bias_mad_pdf_output),
              mgd.OutputFile(params_mad_pdf_output)],
              mgd.InputFile(all_metrics_file),
              config['plot_mad_threshold'],
              config['plot_numreads_threshold'],
              config['plot_quality_threshold']
            )
    )



    return workflow
