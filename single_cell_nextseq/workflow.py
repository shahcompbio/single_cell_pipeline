import os
import pypeliner
import pypeliner.managed as mgd

import single_cell_nextseq.tasks


scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
collect_metrics_script = os.path.join(scripts_directory, 'collect_metrics.py')
extract_quality_metrics_script = os.path.join(scripts_directory, 'extract_quality_metrics.py')


def create_summary_workflow(hmm_segments, hmm_reads, hmm_metrics, metrics_summary, gc_matrix, sample_info_filename, cn_matrix, config, args):#, sample_info_file, sample_id, config):


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
    plot_hmmcopy_script = os.path.join(scripts_directory, 'plot_hmmcopy.py')

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
        func=single_cell_nextseq.tasks.concatenate_csv,
        args=(
            hmm_segments,
            mgd.OutputFile(hmmcopy_segments_filename),
        ),
    )

    workflow.transform(
        name='merge_reads',
        func=single_cell_nextseq.tasks.concatenate_csv,
        args=(
            hmm_reads,
            mgd.OutputFile(hmmcopy_reads_filename),
        ),
    )

    workflow.transform(
        name='merge_hmm_metrics',
        func=single_cell_nextseq.tasks.concatenate_csv,
        args=(
            hmm_metrics,
            mgd.OutputFile(hmmcopy_hmm_metrics_filename),
        ),
    )

    workflow.transform(
        name='merge_summary_metrics',
        func=single_cell_nextseq.tasks.concatenate_csv,
        args=(
            metrics_summary,
            mgd.OutputFile(metrics_summary_filename),
        ),
    )

    workflow.transform(
        name='merge_gc_metrics',
        func=single_cell_nextseq.tasks.merge_csv,
        args=(
            gc_matrix,
            mgd.OutputFile(gc_metrics_filename),
            'outer',
            'gc'
        ),
    )

    workflow.transform(
        name='merge_cn_metrics',
        func=single_cell_nextseq.tasks.merge_csv,
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




def create_fastqc_workflow(fastq_1, fastq_2, trim_1_trim, trim_2_trim, config, metrics_directory, sample_id):#, fastq_1_basename, fastq_2_basename):

    fastqc_1_html_template = os.path.join(metrics_directory, 'fastqc', '{sample_id}_R1.fastqc.html')
    fastqc_1_zip_template = os.path.join(metrics_directory, 'fastqc', '{sample_id}_R1.fastqc.zip')
    fastqc_2_html_template = os.path.join(metrics_directory, 'fastqc', '{sample_id}_R2.fastqc.html')
    fastqc_2_zip_template = os.path.join(metrics_directory, 'fastqc', '{sample_id}_R2.fastqc.zip')

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='produce_fastqc_report_1',
        func=single_cell_nextseq.tasks.produce_fastqc_report,
        args=(
            mgd.InputFile(fastq_1),
            mgd.OutputFile('fastqc_1_html'),
            mgd.OutputFile('fastqc_1_plots'),
            mgd.TempSpace('fastqc_1_temp'),
            config['fastqc']
        ),
    )

    workflow.transform(
        name='produce_fastqc_report_2',
        func=single_cell_nextseq.tasks.produce_fastqc_report,
        args=(
            mgd.InputFile(fastq_2),
            mgd.OutputFile('fastqc_2_html'),
            mgd.OutputFile('fastqc_2_plots'),
            mgd.TempSpace('fastqc_2_temp'),
            config['fastqc']
        ),
    )


    workflow.transform(
        name='run_trimgalore',
        # axes=('sample_id',),
        func=single_cell_nextseq.tasks.run_trimgalore,
        args=(
            fastq_1,
            fastq_2,
            mgd.OutputFile(trim_1_trim),
            mgd.OutputFile(trim_2_trim),
            mgd.TempOutputFile('fastq_fqrep_1'),
            mgd.TempOutputFile('fastq_fqrep_2'),
            mgd.TempOutputFile('fastq_zip_1'),
            mgd.TempOutputFile('fastq_zip_2'),
            mgd.TempOutputFile('fastq_rep_1'),
            mgd.TempOutputFile('fastq_rep_2'),
            mgd.TempSpace('trim_temp'),
            config
        ),
    )

    return workflow

def create_alignment_workflow(
    fastq_1_filename,
    fastq_2_filename,
    bam_filename,
    bam_index_filename,
    ref_genome,
    read_group,
    metrics_summary_filename,
    gc_matrix_filename,
    metrics_directory,
    sample_id,
    samplesheet,
    config):
    gc_metrics_script = os.path.join(scripts_directory, 'gen_cn_matrix.py')


    markdups_metrics_filename = os.path.join(metrics_directory, 'markdups_metrics/{}.markdups_metrics.txt'.format(sample_id))
    flagstat_metrics_filename = os.path.join(metrics_directory, 'flagstat_metrics/{}.flagstat_metrics.txt'.format(sample_id))
    wgs_metrics_filename = os.path.join(metrics_directory, 'wgs_metrics/{}.wgs_metrics.txt'.format(sample_id))
    gc_metrics_filename = os.path.join(metrics_directory, 'gc_metrics/{}.gc_metrics.txt'.format(sample_id))
    gc_summary_filename = os.path.join(metrics_directory, 'gc_metrics/{}.gc_metrics.summ.txt'.format(sample_id))
    gc_chart_filename = os.path.join(metrics_directory, 'gc_metrics/{}.gc_metrics.pdf'.format(sample_id))
    insert_metrics_filename = os.path.join(metrics_directory, 'insert_metrics/{}.insert_metrics.txt'.format(sample_id))
    insert_histogram_filename = os.path.join(metrics_directory, 'insert_metrics/{}.insert_metrics.pdf'.format(sample_id))

    workflow = pypeliner.workflow.Workflow()

    workflow.commandline(
        name='aln_read_1',
        ctx={'mem': 6},
        args=(
            config['bwa'],
            'aln',
            mgd.InputFile(ref_genome),
            mgd.InputFile(fastq_1_filename),
            '>',
            mgd.TempOutputFile('read_1.sai'),
        ),
    )

    workflow.commandline(
        name='aln_read_2',
        ctx={'mem': 6},
        args=(
            config['bwa'], 'aln',
            mgd.InputFile(ref_genome),
            mgd.InputFile(fastq_2_filename),
            '>',
            mgd.TempOutputFile('read_2.sai'),
        ),
    )

    workflow.commandline(
        name='sampe',
        ctx={'mem': 6},
        args=(
            config['bwa'], 'sampe',
            '-r', read_group,
            mgd.InputFile(ref_genome),
            mgd.TempInputFile('read_1.sai'),
            mgd.TempInputFile('read_2.sai'),
            mgd.InputFile(fastq_1_filename),
            mgd.InputFile(fastq_2_filename),
            '|',
            config['samtools'], 'view',
            '-bSh', '-',
            '>',
            mgd.TempOutputFile('aligned.bam'),
        ),
    )

    workflow.transform(
        name='bam_sort',
        ctx={'mem': 16},
        func=single_cell_nextseq.tasks.bam_sort,
        args=(
            mgd.TempInputFile('aligned.bam'),
            mgd.TempOutputFile('sorted.bam'),
            config
        ),
    )

    workflow.transform(
        name='bam_markdups',
        ctx={'mem': 16},
        func=single_cell_nextseq.tasks.bam_markdups,
        args=(
            mgd.TempInputFile('sorted.bam'),
            mgd.OutputFile(bam_filename),
            mgd.OutputFile(markdups_metrics_filename),
            config
        ),
    )

    workflow.commandline(
        name='bam_index',
        ctx={'mem': 16},
        args=(
            config['samtools'], 'index',
            mgd.InputFile(bam_filename),
            # mgd.OutputFile(bam_index_filename),
        ),
    )

    workflow.commandline(
        name='bam_flagstat',
        ctx={'mem': 16},
        args=(
            config['samtools'], 'flagstat',
            mgd.InputFile(bam_filename),
            '>',
            mgd.OutputFile(flagstat_metrics_filename),
        ),
    )

    workflow.transform(
        name='bam_collect_wgs_metrics',
        ctx={'mem': 24},
        func=single_cell_nextseq.tasks.bam_collect_wgs_metrics,
        args=(
            mgd.InputFile(bam_filename),
            mgd.InputFile(ref_genome),
            mgd.OutputFile(wgs_metrics_filename),
            config,
        ),
    )

    workflow.transform(
        name='bam_collect_gc_metrics',
        ctx={'mem': 24},
        func=single_cell_nextseq.tasks.bam_collect_gc_metrics,
        args=(
            mgd.InputFile(bam_filename),
            mgd.InputFile(ref_genome),
            mgd.OutputFile(gc_metrics_filename),
            mgd.OutputFile(gc_summary_filename),
            mgd.OutputFile(gc_chart_filename),
            config
        ),
    )

    workflow.transform(
        name='bam_collect_insert_metrics',
        ctx={'mem': 24},
        func=single_cell_nextseq.tasks.bam_collect_insert_metrics,
        args=(
            mgd.InputFile(bam_filename),
            mgd.InputFile(flagstat_metrics_filename),
            mgd.OutputFile(insert_metrics_filename),
            mgd.OutputFile(insert_histogram_filename),
            config
        ),
    )
    
    workflow.commandline(
        name='collect_metrics',
        ctx={'mem': 16},
        args=(
            config['python'],
            collect_metrics_script,
            mgd.InputFile(flagstat_metrics_filename),
            mgd.InputFile(markdups_metrics_filename),
            mgd.InputFile(insert_metrics_filename),
            mgd.InputFile(wgs_metrics_filename),
            mgd.OutputFile(metrics_summary_filename),
            mgd.InputFile(samplesheet),
            sample_id,
        ),
    )

    workflow.commandline(
        name='collect_gc_metrics',
        ctx={'mem': 16},
        args=(
            config['python'],
            gc_metrics_script,
            '--separator', 'comma',
            '--input', mgd.InputFile(gc_metrics_filename),
            '--output', mgd.OutputFile(gc_matrix_filename),
            '--sample_id', sample_id,
            '--type', 'gcbias', 
            '--column_name', 'NORMALIZED_COVERAGE'
        ),
    )


    return workflow


def create_hmmcopy_workflow(
    bam_filename,
    wig_filename,
    corrected_reads_filename,
    segments_filename,
    parameters_filename,
    posterior_marginals_filename,
    hmm_metrics_filename,
    cnmatrix_filename,
    sample_id,
    config):

    cn_metrics_script = os.path.join(scripts_directory, 'gen_cn_matrix.py')

    chromosomes = config['chromosomes']

    workflow = pypeliner.workflow.Workflow()

    workflow.commandline(
        name='count_reads',
        ctx={'mem': 4},
        args=(
            config['readcounter'],
            '-w', str(config['bin_size']),
            '-q', str(config['min_mqual']),
            '-c', ','.join(chromosomes),
            mgd.InputFile(bam_filename),
            '>',
            mgd.OutputFile(wig_filename),
        ),
    )

    workflow.transform(
        name='run_hmmcopy',
        ctx={'mem': 4},
        func=single_cell_nextseq.tasks.run_hmmcopy,
        args=(
            mgd.InputFile(wig_filename),
            mgd.OutputFile(corrected_reads_filename),
            mgd.OutputFile(segments_filename),
            mgd.OutputFile(parameters_filename),
            mgd.OutputFile(posterior_marginals_filename),
            sample_id,
            config,
        ),
    )

    workflow.commandline(
        name='extract_quality_metrics',
        ctx={'mem': 4},
        args=(
            config['python'],
            extract_quality_metrics_script,
            '--hmmcopy_params', mgd.InputFile(parameters_filename),
            '--hmmcopy_corrected_reads', mgd.InputFile(corrected_reads_filename),
            '--hmmcopy_segments', mgd.InputFile(segments_filename),
            '--out_file', mgd.OutputFile(hmm_metrics_filename),
            '--sample_id', sample_id,
        ),
    )


    workflow.commandline(
        name='collect_gc_metrics',
        ctx={'mem': 16},
        args=(
            config['python'],
            cn_metrics_script,
            '--separator', 'comma',
            '--input', mgd.InputFile(corrected_reads_filename),
            '--output', mgd.OutputFile(cnmatrix_filename),
            '--sample_id', sample_id,
            '--type', 'hmmcopy_corrected_reads', 
            '--column_name', 'integer_copy_number'
        ),
    )

    return workflow

