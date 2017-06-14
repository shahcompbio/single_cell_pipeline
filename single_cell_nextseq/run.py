import os
import sys
import argparse

import pypeliner
import pypeliner.managed as mgd

import single_cell_nextseq.tasks
import single_cell_nextseq.workflow


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
        with open(args.config_file) as file:
            config = yaml.load(file)

    except IOError as e:
        print 'Unable to open config file: {0}'.format(args.config_file)

    sample_sheet_filename = os.path.join(args['nextseq_dir'], 'SampleSheet.csv')

    try:
        with open(sample_sheet_file) as file:
            lines = [x.strip('\n').strip(',') for x in file.readlines()]
        
        run_id = [s.split(',')[1] for s in lines if 'Experiment Name,' in s][0]
        
        library_id = [s.split(',')[1] for s in lines if 'Description,' in s][0]
        
        start_index = lines.index('[Data]')+2

        sample_ids = []
        fastq_1_basename = {}
        fastq_2_basename = {}
        for i, line in zip(range(num_samples), lines[start_index:]):
            sample_id = line.split(',')[0]
            fastq_1_basename[sample_id] = '{0}_S{1}_R1_001'
            fastq_2_basename[sample_id] = '{0}_S{1}_R2_001'
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

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('fastq_1_basename', 'sample_id', axes_origin=[]),
        value=fastq_1_basename,
        axes=('sample_id',),
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('fastq_2_basename', 'sample_id', axes_origin=[]),
        value=fastq_2_basename,
        axes=('sample_id',),
    )

    workflow.transform(
        name='demultiplex_fastq_files',
        func=single_cell_nextseq.tasks.demultiplex_fastq_files,
        args=(
            mgd.InputFile(sample_sheet_filename),
            args['nextseq_dir'],
            mgd.TempOutputFile('fastq_1', 'sample_id', axes_origin=[]),
            mgd.TempOutputFile('fastq_2', 'sample_id', axes_origin=[]),
            mgd.TempSpace('demultiplex_temp'),
        ),
    )

    workflow.transform(
        name='produce_fastqc_report_1',
        axes=('sample_id',),
        func=single_cell_nextseq.tasks.produce_fastqc_report,
        args=(
            mgd.TempInputFile('fastq_1', 'sample_id'),
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
            mgd.TempInputFile('fastq_2', 'sample_id'),
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
            mgd.TempInputFile('fastq_1', 'sample_id'),
            mgd.TempInputFile('fastq_2', 'sample_id'),
            mgd.TempOutputFile('fastq_trim_1', 'sample_id'),
            mgd.TempOutputFile('fastq_trim_2', 'sample_id'),
            mgd.Template(trimgalore_results_template, 'sample_id'),
            mgd.TempInputObj('fastq_1_basename', 'sample_id'),
            mgd.TempInputObj('fastq_2_basename', 'sample_id'),
            config['adapter'],
            config['adapter2'],
        ),
    )

    workflow.subworkflow(
        name='alignment_workflow',
        axes=('sample_id',),
        func=create_alignment_workflow,
        args=(
            mgd.TempInputFile('fastq_trim_1', 'sample_id'),
            mgd.TempInputFile('fastq_trim_2', 'sample_id'),
            mgd.TempOutputFile('bam', 'sample_id'),
            mgd.TempOutputFile('bam_index', 'sample_id'),
            mgd.InputFile(ref_genome),
            mgd.Template(read_group_template, 'sample_id'),
            mgd.OutputFile(metrics_summary_template, 'sample_id'),
            metrics_directory,
            config,
        ),
    )

    pyp.run(workflow)

