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
            'bwa', 'mem', '-C', '-M', '-R', readgroup_literal,
            reference, fastq1, fastq2,
            '>', output,
            **kwargs)
    except pypeliner.commandline.CommandLineException:
        pypeliner.commandline.execute(
            'bwa', 'mem', '-C', '-M', '-R', readgroup,
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
            'bwa', 'sampe', '-r', readgroup_literal, reference, read_1_sai,
            read_2_sai, fastq1, fastq2, '>', output,
            **kwargs)
    except pypeliner.commandline.CommandLineException:
        pypeliner.commandline.execute(
            'bwa', 'sampe', '-r', readgroup, reference, read_1_sai,
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


def biobloom_categorizer(fastq1, fastq2, tempdir, biobloom_count_metrics, disable_biobloom, docker_image,
                         biobloom_filters, ref_type):
    tempdir = os.path.join(tempdir, 'biobloom')
    if not os.path.exists(tempdir):
        helpers.makedirs(tempdir)

    cmd = [
        "biobloomcategorizer",
        "--fq",
        "-e",
        "-p",
        tempdir + "/biobloom",
        "-f",
        " ".join(biobloom_filters),
        fastq1,
        fastq2,
    ]

    if not disable_biobloom:
        pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    extract_biobloom_metrics(tempdir, biobloom_count_metrics, disable_biobloom, biobloom_filters)

    if disable_biobloom:
        final_fastq1 = fastq1
        final_fastq2 = fastq2
    elif ref_type == 'grch37':
        final_fastq1 = os.path.join(tempdir, "biobloom_GRCh37-lite_1.fq")
        final_fastq2 = os.path.join(tempdir, "biobloom_GRCh37-lite_2.fq")
    elif ref_type == 'mm10':
        final_fastq1 = os.path.join(tempdir, "biobloom_mm10_build38_mouse_1.fq")
        final_fastq2 = os.path.join(tempdir, "biobloom_mm10_build38_mouse_2.fq")
    else:
        raise Exception("Unknown reference type specified")

    return final_fastq1, final_fastq2


def get_num_reads(fastq1, fastq2):
    count_1 = 0
    count_2 = 0
    with open(fastq1) as fastqdata:
        for _ in fastqdata:
            count_1 += 1

    with open(fastq2) as fastqdata:
        for _ in fastqdata:
            count_2 += 1

    if not count_1 == count_2:
        ValueError('Two biobloom FastQ counts are not matching')
    else:
        return (count_1 / 4) + (count_2 / 4)


def extract_biobloom_metrics(tempdir, path, disable_biobloom, biobloom_filters):
    biobloom_filters = [os.path.basename(val).replace('.bf', '') for val in biobloom_filters]
    biobloom_filters = ['biobloom_' + val for val in biobloom_filters]

    biobloom_filters.append('biobloom_multiMatch')
    biobloom_filters.append('biobloom_noMatch')

    counts = {}
    for biobloom_filter in biobloom_filters:
        if disable_biobloom:
            counts[biobloom_filter] = 'NA'
        else:
            numreads = get_num_reads(
                os.path.join(tempdir, '{}_1.fq'.format(biobloom_filter)),
                os.path.join(tempdir, '{}_2.fq'.format(biobloom_filter)),
            )
            counts[biobloom_filter] = numreads

    writer = open(path, 'w')
    writer.write(','.join(counts.keys()) + '\n')
    writer.write(','.join(str(v) for v in counts.values()))
    writer.close()


def add_comment_bam_header(infile, outfile, comment):
    with pysam.AlignmentFile(infile, mode='r', check_sq=False) as inbam:
        header = inbam.header.to_dict()
        header['CO'] = comment

        with pysam.AlignmentFile(outfile, header=header, mode='wh') as outbam:
            for read in inbam.fetch(until_eof=True):
                outbam.write(read)
