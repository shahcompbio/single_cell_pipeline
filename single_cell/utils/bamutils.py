'''
Created on Feb 19, 2018

@author: dgrewal
'''
import pypeliner
import shutil
import os
from helpers import makedirs


def produce_fastqc_report(fastq_filename, output_html, output_plots, temp_dir):
    makedirs(temp_dir)

    pypeliner.commandline.execute(
        'fastqc',
        '--outdir=' + temp_dir,
        fastq_filename)


    fastq_basename = os.path.basename(fastq_filename)
    if fastq_basename.endswith(".fastq.gz"):
        fastq_basename = fastq_basename.replace(".fastq.gz", "")
    elif fastq_basename.endswith(".fq.gz"):
        fastq_basename = fastq_basename.replace(".fq.gz", "")
    else:
        raise Exception("Unknown file type")

    output_basename = os.path.join(temp_dir, fastq_basename)

    shutil.move(output_basename + '_fastqc.zip', output_plots)
    shutil.move(output_basename + '_fastqc.html', output_html)


def bwa_mem_paired_end(fastq1, fastq2, output,
                         reference, readgroup):
    """
    run bwa aln on both fastq files,
    bwa sampe to align, and convert to bam with samtools view
    """
    pypeliner.commandline.execute(
        'bwa', 'mem', '-M', '-R', readgroup,
        reference, fastq1, fastq2, '|',
        'samtools', 'view', '-bSh', '-',
        '>', output,
    )


def bwa_aln_paired_end(fastq1, fastq2, output, tempdir,
                         reference, readgroup):
    """
    run bwa aln on both fastq files,
    bwa sampe to align, and convert to bam with samtools view
    """
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    read_1_sai = os.path.join(tempdir, 'read_1.sai')
    read_2_sai = os.path.join(tempdir, 'read_2.sai')

    pypeliner.commandline.execute(
        'bwa',
        'aln',
        reference,
        fastq1,
        '>',
        read_1_sai
    )

    pypeliner.commandline.execute(
        'bwa',
        'aln',
        reference,
        fastq2,
        '>',
        read_2_sai,

    )

    pypeliner.commandline.execute(
        'bwa', 'sampe',
        '-r', readgroup,
        reference,
        read_1_sai,
        read_2_sai,
        fastq1,
        fastq2,
        '|',
        'samtools', 'view',
        '-bSh', '-',
        '>',
        output,
    )

def bam_index(infile, outfile):
    pypeliner.commandline.execute(
        'samtools', 'index',
        infile,
        outfile
        )


def bam_flagstat(bam, metrics):

    pypeliner.commandline.execute(
        'samtools', 'flagstat',
        bam,
        '>',
        metrics
    )


def bam_merge(bams, output, region=None):

    cmd = ['samtools', 'merge']
    if region:
        cmd.extend(['-R', region])

    cmd.append(output)
    cmd.extend(bams)

    pypeliner.commandline.execute(*cmd)

