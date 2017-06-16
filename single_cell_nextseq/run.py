import os
import sys
import yaml
import argparse

import pypeliner
import pypeliner.managed as mgd

import single_cell_nextseq.tasks
import single_cell_nextseq.workflow


scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
plot_hmmcopy_script = os.path.join(scripts_directory, 'plot_hmmcopy.py')


if __name__ == '__main__':
    main()


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    pypeliner.app.add_arguments(parser)

    parser.add_argument('nextseq_dir',
                        help='''Path to input nextseq directory.''')

    parser.add_argument('out_dir',
                        help='''Path to output files.''')

    parser.add_argument('config_file',
                        help='''Path to yaml config file.''')

    args = vars(parser.parse_args())
    args['tmpdir'] = os.path.join(args['out_dir'], 'tmp')

    pyp = pypeliner.app.Pypeline(config=args)

    try:
        with open(args['config_file']) as file:
            config = yaml.load(file)

    except IOError as e:
        print 'Unable to open config file: {0}'.format(args['config_file'])

    sample_sheet_filename = os.path.join(args['nextseq_dir'], 'SampleSheet.csv')
    fastq_directory = os.path.join(args['out_dir'], 'fastq')

    try:
        with open(sample_sheet_filename) as file:
            lines = [x.strip('\n').strip(',') for x in file.readlines()]
        
        run_id = [s.split(',')[1] for s in lines if 'Experiment Name,' in s][0]
        
        library_id = [s.split(',')[1] for s in lines if 'Description,' in s][0]
        
        start_index = lines.index('[Data]')+2

        num_samples = len(lines[start_index:])

        sample_ids = []
        fastq_1_basenames = {}
        fastq_2_basenames = {}
        fastq_1_filenames = {}
        fastq_2_filenames = {}
        for i, line in zip(range(num_samples), lines[start_index:]):
            sample_id = line.split(',')[0]
            fastq_1_basenames[sample_id] = '{0}_S{1}_R1_001'.format(sample_id, str(i+1))
            fastq_2_basenames[sample_id] = '{0}_S{1}_R2_001'.format(sample_id, str(i+1))
            fastq_1_filenames[sample_id] = os.path.join(fastq_directory, fastq_1_basenames[sample_id] + '.fastq.gz')
            fastq_2_filenames[sample_id] = os.path.join(fastq_directory, fastq_2_basenames[sample_id] + '.fastq.gz')
            sample_ids.append(sample_id)

    except IOError as e:
        print 'Unable to open file \'SampleSheet.csv\' in directory: {0}'.format(args['nextseq_dir'])

    if 'read_group' in config.keys():
        read_group_template = (
            '@RG\tID:' + str(library_id) + '_{sample_id}_' + str(run_id) + 
            '\tPL:' + str(config['read_group']['PL']) +
            '\tPU:' + str(run_id) +
            '\tLB:' + str(library_id) + '_{sample_id}' +
            '\tSM:' + '{sample_id}' +
            '\tCN:' + str(config['read_group']['CN']))
    
    else:
        warnings.warn('Config file does not contain read group information! ' + 
                      'This will affect duplicate marking if BAMs are later merged. ' +
                      'Creating BAM without read group information in header.')

    trimgalore_results_template = os.path.join(args['out_dir'], 'trim', '{sample_id}')
    metrics_directory = os.path.join(args['out_dir'], 'metrics')
    fastqc_1_html_template = os.path.join(metrics_directory, 'fastqc', '{sample_id}_R1.fastqc.html')
    fastqc_1_zip_template = os.path.join(metrics_directory, 'fastqc', '{sample_id}_R1.fastqc.zip')
    fastqc_2_html_template = os.path.join(metrics_directory, 'fastqc', '{sample_id}_R2.fastqc.html')
    fastqc_2_zip_template = os.path.join(metrics_directory, 'fastqc', '{sample_id}_R2.fastqc.zip')
    metrics_summary_template = os.path.join(metrics_directory, 'summary', '{sample_id}_summary.csv')
    metrics_summary_filename = os.path.join(metrics_directory, 'summary', 'summary.csv')

    bam_directory = os.path.join(args['out_dir'], 'bams')
    bam_template = os.path.join(bam_directory, '{sample_id}.bam')
    bam_index_template = os.path.join(bam_directory, '{sample_id}.bam.bai')

    hmmcopy_directory = os.path.join(args['out_dir'], 'hmmcopy')
    hmmcopy_wig_template = os.path.join(hmmcopy_directory, 'hmmcopy', '{sample_id}_readcount.wig')
    hmmcopy_reads_template = os.path.join(hmmcopy_directory, 'hmmcopy', '{sample_id}_reads.csv')
    hmmcopy_segments_template = os.path.join(hmmcopy_directory, 'hmmcopy', '{sample_id}_segments.csv')
    hmmcopy_parameters_template = os.path.join(hmmcopy_directory, 'hmmcopy', '{sample_id}_parameters.csv')
    hmmcopy_posteriors_template = os.path.join(hmmcopy_directory, 'hmmcopy', '{sample_id}_posteriors.csv')
    hmmcopy_hmm_metrics_template = os.path.join(hmmcopy_directory, 'hmmcopy', '{sample_id}_hmm_metrics.csv')

    hmmcopy_segments_filename = os.path.join(hmmcopy_directory, 'hmmcopy', 'segments.csv')
    hmmcopy_reads_filename = os.path.join(hmmcopy_directory, 'hmmcopy', 'reads.csv')
    hmmcopy_hmm_metrics_filename = os.path.join(hmmcopy_directory, 'hmmcopy', 'hmm_metrics.csv')

    plots_directory = os.path.join(args['out_dir'], 'plots')
    reads_plot_filename = os.path.join(plots_directory, 'corrected_reads.pdf')
    bias_plot_filename = os.path.join(plots_directory, 'bias.pdf')
    segs_plot_filename = os.path.join(plots_directory, 'segments.pdf')

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('fastq_1_basename', 'sample_id', axes_origin=[]),
        value=fastq_1_basenames,
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('fastq_2_basename', 'sample_id', axes_origin=[]),
        value=fastq_2_basenames,
    )

    workflow.transform(
        name='demultiplex_fastq_files',
        ctx={'mem': 64},
        func=single_cell_nextseq.tasks.demultiplex_fastq_files,
        args=(
            mgd.InputFile(sample_sheet_filename),
            args['nextseq_dir'],
            mgd.OutputFile('fastq_1', 'sample_id', fnames=fastq_1_filenames, axes_origin=[]),
            mgd.OutputFile('fastq_2', 'sample_id', fnames=fastq_2_filenames, axes_origin=[]),
            mgd.TempSpace('demultiplex_temp'),
        ),
    )

    workflow.transform(
        name='produce_fastqc_report_1',
        axes=('sample_id',),
        func=single_cell_nextseq.tasks.produce_fastqc_report,
        args=(
            mgd.InputFile('fastq_1', 'sample_id', fnames=fastq_1_filenames),
            mgd.TempInputObj('fastq_1_basename', 'sample_id'),
            mgd.OutputFile('fastqc_1_html', 'sample_id', template=fastqc_1_html_template),
            mgd.OutputFile('fastqc_1_plots', 'sample_id', template=fastqc_1_zip_template),
            mgd.TempSpace('fastqc_1_temp', 'sample_id'),
        ),
    )

    workflow.transform(
        name='produce_fastqc_report_2',
        axes=('sample_id',),
        func=single_cell_nextseq.tasks.produce_fastqc_report,
        args=(
            mgd.InputFile('fastq_2', 'sample_id', fnames=fastq_2_filenames),
            mgd.TempInputObj('fastq_2_basename', 'sample_id'),
            mgd.OutputFile('fastqc_2_html', 'sample_id', template=fastqc_2_html_template),
            mgd.OutputFile('fastqc_2_plots', 'sample_id', template=fastqc_2_zip_template),
            mgd.TempSpace('fastqc_2_temp', 'sample_id'),
        ),
    )

    workflow.transform(
        name='run_trimgalore',
        axes=('sample_id',),
        func=single_cell_nextseq.tasks.run_trimgalore,
        args=(
            mgd.InputFile('fastq_1', 'sample_id', fnames=fastq_1_filenames),
            mgd.InputFile('fastq_2', 'sample_id', fnames=fastq_2_filenames),
            mgd.TempInputObj('fastq_1_basename', 'sample_id'),
            mgd.TempInputObj('fastq_2_basename', 'sample_id'),
            mgd.TempOutputFile('fastq_trim_1', 'sample_id'),
            mgd.TempOutputFile('fastq_trim_2', 'sample_id'),
            mgd.Template(trimgalore_results_template, 'sample_id'),
            config['adapter'],
            config['adapter2'],
        ),
    )

    workflow.subworkflow(
        name='alignment_workflow',
        axes=('sample_id',),
        func=single_cell_nextseq.workflow.create_alignment_workflow,
        args=(
            mgd.TempInputFile('fastq_trim_1', 'sample_id'),
            mgd.TempInputFile('fastq_trim_2', 'sample_id'),
            mgd.OutputFile('bam', 'sample_id', template=bam_template),
            mgd.OutputFile('bam_index', 'sample_id', template=bam_index_template),
            mgd.InputFile(config['ref_genome']),
            mgd.Template(read_group_template, 'sample_id'),
            mgd.OutputFile(metrics_summary_template, 'sample_id'),
            metrics_directory,
            mgd.InputInstance('sample_id'),
            config,
        ),
    )

    workflow.subworkflow(
        name='hmmcopy_workflow',
        axes=('sample_id',),
        func=single_cell_nextseq.workflow.create_hmmcopy_workflow,
        args=(
            mgd.InputFile('bam', 'sample_id', template=bam_template),
            mgd.OutputFile(hmmcopy_wig_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_reads_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_segments_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_parameters_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_posteriors_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_hmm_metrics_template, 'sample_id'),
            mgd.InputInstance('sample_id'),
            config,
        ),
    )

    workflow.transform(
        name='merge_tables',
        func=single_cell_nextseq.tasks.concatenate_csv,
        args=(
            mgd.InputFile(hmmcopy_segments_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_segments_filename),
        ),
    )

    workflow.transform(
        name='merge_reads',
        func=single_cell_nextseq.tasks.concatenate_csv,
        args=(
            mgd.InputFile(hmmcopy_reads_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_reads_filename),
        ),
    )

    workflow.transform(
        name='merge_hmm_metrics',
        func=single_cell_nextseq.tasks.concatenate_csv,
        args=(
            mgd.InputFile(hmmcopy_hmm_metrics_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_hmm_metrics_filename),
        ),
    )

    workflow.transform(
        name='merge_summary_metrics',
        func=single_cell_nextseq.tasks.concatenate_csv,
        args=(
            mgd.InputFile(metrics_summary_template, 'sample_id'),
            mgd.OutputFile(metrics_summary_filename),
        ),
    )

    workflow.commandline(
        name='plot_hmm_copy',
        args=(
            'python',
            plot_hmmcopy_script,
            '--corrected_reads', mgd.InputFile(hmmcopy_reads_filename),
            '--segments', mgd.InputFile(hmmcopy_segments_filename),
            '--hmm_metrics', mgd.InputFile(hmmcopy_hmm_metrics_filename),
            '--align_metrics', mgd.InputFile(metrics_summary_filename),
            '--ref_genome', mgd.InputFile(config['ref_genome']),
            '--num_states', config['num_states'],
            '--reads_output', mgd.OutputFile(reads_plot_filename),
            '--bias_output', mgd.OutputFile(bias_plot_filename),
            '--segs_output', mgd.OutputFile(segs_plot_filename),
            '--plot_title', 'QC pipeline metrics',
        ),
    )

    pyp.run(workflow)

