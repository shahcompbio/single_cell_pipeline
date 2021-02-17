'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import shutil

import pypeliner
import pysam
from single_cell.utils import helpers

from single_cell.utils.helpers import makedirs


def produce_fastqc_report(fastq_filename, output_html, output_plots, temp_dir,):
    makedirs(temp_dir)

    pypeliner.commandline.execute(
        'fastqc',
        '--outdir=' + temp_dir,
        fastq_filename,
        )

    fastq_basename = os.path.basename(fastq_filename)
    if fastq_basename.endswith(".fastq.gz"):
        fastq_basename = fastq_basename.replace(".fastq.gz", "")
    elif fastq_basename.endswith(".fq.gz"):
        fastq_basename = fastq_basename.replace(".fq.gz", "")
    elif fastq_basename.endswith(".fq"):
        fastq_basename = fastq_basename.replace(".fq", "")
    elif fastq_basename.endswith(".fastq"):
        fastq_basename = fastq_basename.replace(".fastq", "")
    else:
        raise Exception("Unknown file type")

    output_basename = os.path.join(temp_dir, fastq_basename)

    shutil.move(output_basename + '_fastqc.zip', output_plots)
    shutil.move(output_basename + '_fastqc.html', output_html)


def bwa_mem_paired_end(fastq1, fastq2, output,
                       reference, readgroup,
                       ):
    """
    run bwa aln on both fastq files,
    bwa sampe to align, and convert to bam with samtools view
    """

    try:
        readgroup_literal = '"' + readgroup + '"'
        pypeliner.commandline.execute(
            'bwa', 'mem', '-C', '-M', '-R', readgroup_literal,
            reference, fastq1, fastq2,
            '>', output,
            )
    except pypeliner.commandline.CommandLineException:
        pypeliner.commandline.execute(
            'bwa', 'mem', '-C', '-M', '-R', readgroup,
            reference, fastq1, fastq2,
            '>', output,
            )


def samtools_sam_to_bam(samfile, bamfile,
                        ):
    pypeliner.commandline.execute(
        'samtools', 'view', '-bSh', samfile,
        '>', bamfile,
        )


def bwa_aln_paired_end(fastq1, fastq2, output, tempdir,
                       reference, readgroup,
                       ):
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
        read_1_sai,
        )

    pypeliner.commandline.execute(
        'bwa',
        'aln',
        reference,
        fastq2,
        '>',
        read_2_sai,
        )

    try:
        readgroup_literal = '"' + readgroup + '"'
        pypeliner.commandline.execute(
            'bwa', 'sampe', '-r', readgroup_literal, reference, read_1_sai,
            read_2_sai, fastq1, fastq2, '>', output,
            )
    except pypeliner.commandline.CommandLineException:
        pypeliner.commandline.execute(
            'bwa', 'sampe', '-r', readgroup, reference, read_1_sai,
            read_2_sai, fastq1, fastq2, '>', output,
            )


def bam_index(infile, outfile):
    pypeliner.commandline.execute(
        'samtools', 'index',
        infile,
        outfile,
        )


def bam_flagstat(bam, metrics):
    pypeliner.commandline.execute(
        'samtools', 'flagstat',
        bam,
        '>',
        metrics,
        )


def bam_merge(bams, output, **kwargs):
    if isinstance(bams, dict):
        bams = bams.values()

    cmd = ['samtools', 'merge', '-f']
    if kwargs.get('region'):
        cmd.extend(['-R', kwargs.get('region')])

    cmd.append(output)
    cmd.extend(bams)

    pypeliner.commandline.execute(*cmd)


def bam_view(bam, output, region):
    cmd = ['samtools', 'view', '-b', bam, '-o', output, region]

    pypeliner.commandline.execute(*cmd)


def add_comment_bam_header(infile, outfile, comment):
    with pysam.AlignmentFile(infile, mode='r', check_sq=False) as inbam:
        header = inbam.header.to_dict()
        header['CO'] = comment

        with pysam.AlignmentFile(outfile, header=header, mode='wh') as outbam:
            for read in inbam.fetch(until_eof=True):
                outbam.write(read)
