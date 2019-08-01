import os

import pypeliner
import pysam
import single_cell.utils.helpers as helpers
import yaml
from single_cell.workflows.lumpy import generate_histogram


def process_bam(
        input_bam, discordant_bam, split_bam, histogram,
        tempdir, docker_image=None, tag=None,
        N=10000, skip=100000, min_elements=1000, mads=10,
        X=4, read_length=101
):
    helpers.makedirs(tempdir)

    discordants = os.path.join(tempdir, "discordants.bam")
    run_samtools_view(input_bam, discordants, docker_image=docker_image)
    if tag:
        sorted_discordants = os.path.join(tempdir, 'discordants.sorted.bam')
        run_samtools_sort(discordants, sorted_discordants, docker_image=docker_image)
        tag_reads(sorted_discordants, discordant_bam, tag)
    else:
        run_samtools_sort(discordants, discordant_bam, docker_image=docker_image)

    split_reads = os.path.join(tempdir, "split_reads.bam")
    run_lumpy_extract_split_reads_bwamem(input_bam, split_reads, docker_image=docker_image)
    if tag:
        sorted_split_reads = os.path.join(tempdir, 'split.sorted.bam')
        run_samtools_sort(split_reads, sorted_split_reads, docker_image=docker_image)
        tag_reads(sorted_split_reads, split_bam, tag)
    else:
        run_samtools_sort(split_reads, split_bam, docker_image=docker_image)

    generate_histogram.gen_histogram(
        input_bam, histogram, N=N, skip=skip, min_elements=min_elements,
        mads=mads, X=X, read_length=read_length
    )


def tag_reads(infile, outfile, sample_id):
    infile = pysam.AlignmentFile(infile, 'rb')
    taggedreads = pysam.AlignmentFile(outfile, "wb", template=infile)
    for read in infile.fetch():
        read.query_name = "{}:{}".format(sample_id, read.query_name)
        taggedreads.write(read)
    infile.close()
    taggedreads.close()


def run_samtools_view(infile, outfile, docker_image=None):
    cmd = ['samtools', 'view', '-b', '-F', '1294', infile, '>', outfile]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def run_lumpy_extract_split_reads_bwamem(infile, outfile, docker_image=None):
    cmd = [
        'samtools', 'view', '-h', infile, '|',
        'extractSplitReads_BwaMem', '-i', 'stdin', '|',
        'samtools', 'view', '-Sb', '-',
        '>', outfile
    ]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def run_samtools_sort(infile, outfile, docker_image=None):
    cmd = ['samtools', 'sort', infile, '-o', outfile]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    cmd = ['samtools', 'index', outfile]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def merge_bams(inputs, output, tempdir, docker_image=None):
    helpers.makedirs(tempdir)

    inputs = inputs.values()

    cmd = ['samtools', 'merge', '-f', output]
    cmd.extend(inputs)

    shell_script_path = os.path.join(tempdir, "run_merge.sh")

    with open(shell_script_path, 'w') as scriptfile:
        scriptfile.write("#!/bin/bash\n")
        cmd = ' '.join(cmd) + '\n'
        scriptfile.write(cmd)

    cmd = ['sh', shell_script_path]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def load_metadata(infile):
    with open(infile) as fileinput:
        data = yaml.safe_load(fileinput)
        return data['mean'], data['stdev']


def run_lumpy(
        tumour_disc, tumour_split, tumour_hist, tumour_mean_stdev, tumour_id,
        normal_disc, normal_split, normal_hist, normal_mean_stdev, normal_id,
        vcf, tempdir, docker_image=None
):
    tumour_mean, tumour_stdev = load_metadata(tumour_mean_stdev)
    normal_mean, normal_stdev = load_metadata(normal_mean_stdev)

    helpers.makedirs(tempdir)
    tempdir = tempdir + '/lumpy'

    tumour_pe = 'id:{},bam_file:{},histo_file:{},mean:{},' \
                'stdev:{},read_length:101,min_non_overlap:101,' \
                'discordant_z:5,back_distance:10,weight:1,' \
                'min_mapping_threshold:20'.format(tumour_id, tumour_disc, tumour_hist, tumour_mean, tumour_stdev)
    tumour_sr = 'id:{},bam_file:{},back_distance:10,weight:1,' \
                'min_mapping_threshold:20'.format(tumour_id, tumour_split)

    normal_pe = 'id:{},bam_file:{},histo_file:{},mean:{},' \
                'stdev:{},read_length:101,min_non_overlap:101,' \
                'discordant_z:5,back_distance:10,weight:1,' \
                'min_mapping_threshold:20'.format(normal_id, normal_disc, normal_hist, normal_mean, normal_stdev)
    normal_sr = 'id:{},bam_file:{},back_distance:10,weight:1,' \
                'min_mapping_threshold:20'.format(normal_id, normal_split)

    cmd = ['lumpy', '-e', '-b', '-mw', 4, '-tt', 0,
           '-pe', tumour_pe, '-sr', tumour_sr,
           '-pe', normal_pe, '-sr', normal_sr,
           '-t', tempdir, '>', vcf]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)
