'''
Created on Jul 6, 2017

@author: dgrewal
'''

import os
import pypeliner
import pypeliner.managed as mgd
import tasks

def create_fastqc_workflow(fastq_r1, fastq_r2, trim_r1, trim_r2, config, lane, sample_id, args, trim):

    scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')

    fastqc_dir = os.path.join(args['out_dir'], 'fastqc')
    html_r1 = os.path.join(fastqc_dir, lane, '{}_R1.html'.format(sample_id))
    html_r2 = os.path.join(fastqc_dir, lane, '{}_R2.html'.format(sample_id))

    plots_r1 = os.path.join(fastqc_dir, lane, '{}_R1.zip'.format(sample_id))
    plots_r2 = os.path.join(fastqc_dir, lane, '{}_R2.zip'.format(sample_id))

    trimgalore_dir = os.path.join(args['out_dir'], 'fastqc', 'trimgalore')
    qc_report_r1 = os.path.join(trimgalore_dir, lane, '{}_QC_R1.html'.format(sample_id))
    qc_report_r2 = os.path.join(trimgalore_dir, lane, '{}_QC_R2.html'.format(sample_id))
    zip_r1 = os.path.join(trimgalore_dir, lane, '{}_R1.zip'.format(sample_id))
    zip_r2 = os.path.join(trimgalore_dir, lane, '{}_R2.zip'.format(sample_id))
    report_r1 = os.path.join(trimgalore_dir, lane, '{}_R1.html'.format(sample_id))
    report_r2 = os.path.join(trimgalore_dir, lane, '{}_R2.html'.format(sample_id))
    

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='produce_fastqc_report_1',
        ctx={'mem': config['low_mem']},
        func=tasks.produce_fastqc_report,
        args=(
            mgd.InputFile(fastq_r1),
            mgd.OutputFile(html_r1),
            mgd.OutputFile(plots_r1),
            mgd.TempSpace('fastqc_1_temp'),
            config['fastqc']
        ),
    )

    workflow.transform(
        name='produce_fastqc_report_2',
        ctx={'mem': config['low_mem']},
        func=tasks.produce_fastqc_report,
        args=(
            mgd.InputFile(fastq_r2),
            mgd.OutputFile(html_r2),
            mgd.OutputFile(plots_r2),
            mgd.TempSpace('fastqc_2_temp'),
            config['fastqc']
        ),
    )

    if trim:
        workflow.commandline(
            name='run_trimgalore',
            ctx={'mem': config['low_mem']},
            args=(
              config['python'],
              os.path.join(scripts_directory, 'run_trimgalore.py'),
                mgd.InputFile(fastq_r1),
                mgd.InputFile(fastq_r2),
                mgd.OutputFile(trim_r1),
                mgd.OutputFile(trim_r2),
                mgd.TempOutputFile(qc_report_r1),
                mgd.TempOutputFile(qc_report_r2),
                mgd.TempOutputFile(zip_r1),
                mgd.TempOutputFile(zip_r2),
                mgd.TempOutputFile(report_r1),
                mgd.TempOutputFile(report_r2),
                mgd.TempSpace('trim_temp'),
                '--adapter', config['adapter'],
                '--adapter2', config['adapter2'],
                '--trimgalore_path', config['trimgalore'],
                '--cutadapt_path', config['cutadapt'],
              )
            )
    else:
        workflow.transform(
            name='copy_files',
            func=tasks.copy_files,
            args=(
                mgd.InputFile(fastq_r1),
                mgd.InputFile(fastq_r2),
                mgd.OutputFile(trim_r1),
                mgd.OutputFile(trim_r2),
            ),
    )

    return workflow

