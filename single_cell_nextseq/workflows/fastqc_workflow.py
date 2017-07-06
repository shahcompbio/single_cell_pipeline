'''
Created on Jul 6, 2017

@author: dgrewal
'''

import os
import pypeliner
import pypeliner.managed as mgd
import single_cell_nextseq.tasks

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

