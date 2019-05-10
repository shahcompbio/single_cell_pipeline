'''
Created on Feb 19, 2018

@author: dgrewal
'''
import pypeliner
import shutil
import os
from helpers import makedirs


def produce_fastqc_report(fastq_filename, output_html, output_plots, temp_dir,
                          **kwargs):
    makedirs(temp_dir)

    pypeliner.commandline.execute(
        'fastqc',
        '--outdir=' + temp_dir,
        fastq_filename,
        **kwargs)

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
                         **kwargs):
    """
    run bwa aln on both fastq files,
    bwa sampe to align, and convert to bam with samtools view
    """

    try:
        readgroup_literal = '"' + readgroup + '"'
        pypeliner.commandline.execute(
            'bwa', 'mem', '-M', '-R', readgroup_literal,
            reference, fastq1, fastq2,
            '>', output,
            **kwargs)
    except pypeliner.commandline.CommandLineException:
        pypeliner.commandline.execute(
            'bwa', 'mem', '-M', '-R', readgroup,
            reference, fastq1, fastq2,
            '>', output,
            **kwargs)


def samtools_sam_to_bam(samfile, bamfile,
                         **kwargs):

    pypeliner.commandline.execute(
        'samtools', 'view', '-bSh', samfile,
        '>', bamfile,
        **kwargs)


def bwa_aln_paired_end(fastq1, fastq2, output, tempdir,
                         reference, readgroup,
                         **kwargs):
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
        **kwargs)

    pypeliner.commandline.execute(
        'bwa',
        'aln',
        reference,
        fastq2,
        '>',
        read_2_sai,
        **kwargs)

    try:
        readgroup_literal = '"' + readgroup + '"'
        pypeliner.commandline.execute(
            'bwa', 'sampe', '-r', readgroup_literal,  reference,  read_1_sai,
            read_2_sai, fastq1, fastq2, '>', output,
            **kwargs)
    except pypeliner.commandline.CommandLineException:
        pypeliner.commandline.execute(
            'bwa', 'sampe', '-r', readgroup,  reference,  read_1_sai,
            read_2_sai, fastq1, fastq2, '>', output,
            **kwargs)


def bam_index(infile, outfile, **kwargs):

    pypeliner.commandline.execute(
        'samtools', 'index',
        infile,
        outfile,
        **kwargs)


def bam_flagstat(bam, metrics, **kwargs):

    pypeliner.commandline.execute(
        'samtools', 'flagstat',
        bam,
        '>',
        metrics,
        **kwargs)


def bam_merge(bams, output, **kwargs):

    if isinstance(bams, dict):
        bams = bams.values()

    cmd = ['samtools', 'merge', '-f']
    if kwargs.get('region'):
        cmd.extend(['-R', kwargs.get('region')])

    cmd.append(output)
    cmd.extend(bams)

    pypeliner.commandline.execute(*cmd, docker_image=kwargs.get('docker_image'))


def bam_view(bam, output, region, **kwargs):

    cmd = ['samtools', 'view', '-b', bam, '-o', output, region]

    pypeliner.commandline.execute(*cmd, **kwargs)


def run_biobloom(fastq1, fastq2, docker_image):

    print "\n\n\n\n\n\n\n\n\n\n\n\nBIOBLOOM\n\n\n\n\n\n\n\n\n\n\n\n"
    print fastq1
    print fastq2
    biobloomcategorizer = "/home/linuxbrew/.linuxbrew/Cellar/biobloomtools/2.3.2/bin/biobloomcategorizer"

    cmd = [
        "biobloomcategorizer",
        "--fq",
        "-e",
        "-p",
        "/biobloom_output/biobloom",
        "-f",
        "/refdata/GCF_002021735.1_Okis_V1_genomic.bf /refdata/GRCh37-lite.bf /refdata/mm10_build38_mouse.bf",
        fastq1,
        fastq2,
    ]

    print cmd

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)
    print "\n\n\n\n\nBIOBLOOM GENERATED\n\n\n\n\n"
    return "./biobloom_output/biobloom_GRCh37-lite_1.fq", "./biobloom_output/biobloom_GRCh37-lite_2.fq"
