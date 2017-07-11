'''
Created on Jul 6, 2017

@author: dgrewal
'''

import os
import pypeliner
import pypeliner.managed as mgd
import tasks

def create_fastqc_workflow(fastq_1, fastq_2, trim_1_trim, trim_2_trim, config, lane, sample_id, args, trim):

    scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')

    hmmcopy_directory = os.path.join(args['out_dir'], 'fastqc')
    html_r1 = os.path.join(hmmcopy_directory, lane, sample_id+'_R1.html')
    html_r2 = os.path.join(hmmcopy_directory, lane, sample_id+'_R2.html')

    plots_r1 = os.path.join(hmmcopy_directory, lane, sample_id+'_R1.zip')
    plots_r2 = os.path.join(hmmcopy_directory, lane, sample_id+'_R2.zip')


    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='produce_fastqc_report_1',
        func=tasks.produce_fastqc_report,
        args=(
            mgd.InputFile(fastq_1),
            mgd.OutputFile(html_r1),
            mgd.OutputFile(plots_r1),
            mgd.TempSpace('fastqc_1_temp'),
            config['fastqc']
        ),
    )

    workflow.transform(
        name='produce_fastqc_report_2',
        func=tasks.produce_fastqc_report,
        args=(
            mgd.InputFile(fastq_2),
            mgd.OutputFile(html_r2),
            mgd.OutputFile(plots_r2),
            mgd.TempSpace('fastqc_2_temp'),
            config['fastqc']
        ),
    )

#     workflow.commandline(
#         name='run_trimgalore',
#         args=(
#           config['python'],
#           os.path.join(scripts_directory, 'run_trimgalore.py'),
#             mgd.InputFile(fastq_1),
#             mgd.InputFile(fastq_2),
#             mgd.OutputFile(trim_1_trim),
#             mgd.OutputFile(trim_2_trim),
#             mgd.TempOutputFile('fastq_fqrep_1'),
#             mgd.TempOutputFile('fastq_fqrep_2'),
#             mgd.TempOutputFile('fastq_zip_1'),
#             mgd.TempOutputFile('fastq_zip_2'),
#             mgd.TempOutputFile('fastq_rep_1'),
#             mgd.TempOutputFile('fastq_rep_2'),
#             mgd.TempSpace('trim_temp'),
#             '--adapter', config['adapter'],
#             '--adapter2', config['adapter2'],
#             '--trimgalore_path', config['trimgalore'],
#             '--cutadapt_path', config['cutadapt'],
#           )
#         )

    if trim:
        workflow.commandline(
            name='run_trimgalore',
            args=(
              config['python'],
              os.path.join(scripts_directory, 'run_trimgalore.py'),
                mgd.InputFile(fastq_1),
                mgd.InputFile(fastq_2),
                mgd.OutputFile(trim_1_trim),
                mgd.OutputFile(trim_2_trim),
                mgd.TempOutputFile('fastq_fqrep_1'),
                mgd.TempOutputFile('fastq_fqrep_2'),
                mgd.TempOutputFile('fastq_zip_1'),
                mgd.TempOutputFile('fastq_zip_2'),
                mgd.TempOutputFile('fastq_rep_1'),
                mgd.TempOutputFile('fastq_rep_2'),
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
                mgd.InputFile(fastq_1),
                mgd.InputFile(fastq_2),
                mgd.OutputFile(trim_1_trim),
                mgd.OutputFile(trim_2_trim),
            ),
         
     
    )

    return workflow

