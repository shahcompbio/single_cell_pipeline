'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import errno
import shutil
import pypeliner

from scripts import RunTrimGalore


def copy_files(inputs,outputs):
    for inp, outp in zip(inputs, outputs):
        shutil.copy(inp, outp)


def makedirs(directory):
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def produce_fastqc_report(fastq_filename, output_html, output_plots, temp_directory, fastqc):
    makedirs(temp_directory)

    pypeliner.commandline.execute(
        fastqc,
        '--outdir=' + temp_directory,
        fastq_filename)

    fastq_basename = os.path.basename(fastq_filename).split('.')[0]
    output_basename = os.path.join(temp_directory, fastq_basename)

    os.rename(output_basename + '_fastqc.zip', output_plots)
    os.rename(output_basename + '_fastqc.html', output_html)
    
def run_trimgalore(seq1, seq2, fq_r1, fq_r2, trimgalore, cutadapt, tempdir,
                   adapter, adapter2, report_r1, report_r2, qc_report_r1,
                   qc_report_r2, qc_zip_r1, qc_zip_r2):

    run_tg = RunTrimGalore(seq1, seq2, fq_r1, fq_r2, trimgalore, cutadapt,
                           tempdir, adapter, adapter2, report_r1, report_r2,
                           qc_report_r1, qc_report_r2, qc_zip_r1, qc_zip_r2)
    run_tg.run_trimgalore()
    run_tg.gather_outputs()
