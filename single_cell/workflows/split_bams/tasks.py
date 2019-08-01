'''
Created on Nov 21, 2017

@author: dgrewal
'''
import gzip
import os

from single_cell.utils import bamutils
from single_cell.utils import helpers

import pypeliner


def split_bam_file_one_job(bam, outbam, regions, samtools_docker, tempdir, ncores=None):
    commands = []
    for region in regions:
        output = outbam[region]
        region = '{}:{}-{}'.format(*region.split('-'))
        commands.append(['samtools', 'view', '-b', bam, '-o', output, region])

    helpers.run_in_gnu_parallel(commands, tempdir, samtools_docker, ncores=ncores)

    commands = []
    for region in regions:
        output = outbam[region]
        commands.append(['samtools', 'index', output, output + ".bai"])

    helpers.run_in_gnu_parallel(commands, tempdir, samtools_docker, ncores=ncores)


def split_bam_file(bam, outbam, interval, samtools_docker):
    outbai = outbam + '.bai'
    bamutils.bam_view(bam, outbam, interval, docker_image=samtools_docker)

    bamutils.bam_index(outbam, outbai, docker_image=samtools_docker)


def split_bam_file_by_reads(bam, outbams, tempspace, intervals, samtools_docker):
    # sort bam by reads and convert to sam

    helpers.makedirs(tempspace)

    headerfile = os.path.join(tempspace, "bam_header.sam")

    cmd = ['samtools', 'view', '-H', bam, '-o', headerfile]
    pypeliner.commandline.execute(*cmd, docker_image=samtools_docker)

    collate_prefix = os.path.join(
        tempspace, os.path.basename(bam) + "_collate_temp"
    )
    collated_bam = os.path.join(tempspace, "bam_file_collated_sam_format.sam")

    cmd = [
        'samtools', 'collate', '-u', '-O', bam, collate_prefix, '|',
        'samtools', 'view', '-', '-o', collated_bam
    ]

    pypeliner.commandline.execute(*cmd, docker_image=samtools_docker)

    tempoutputs = [
        os.path.join(tempspace, os.path.basename(outbams[interval]) + ".split.temp")
        for interval in intervals
    ]

    split(collated_bam, tempoutputs, headerfile=headerfile)

    for inputsam, interval in zip(tempoutputs, intervals):
        outputbam = outbams[interval]

        cmd = ['samtools', 'view', '-Sb', inputsam, '-o', outputbam]

        pypeliner.commandline.execute(*cmd, docker_image=samtools_docker)


def get_file_handle(filename, mode="r"):
    if os.path.splitext(filename)[-1] == ".gz":
        return gzip.open(filename, mode + "b")
    else:
        return open(filename, mode)


def get_file_linecount(infile):
    def blocks(files, size=65536):
        while True:
            b = files.read(size)
            if not b:
                break
            yield b

    with get_file_handle(infile, mode="r") as inputdata:
        numlines = sum(bl.count("\n") for bl in blocks(inputdata))

    return numlines


def split(infile, outfiles, headerfile=None):
    if headerfile:
        with open(headerfile) as headerdata:
            header = headerdata.readlines()

    with get_file_handle(infile, mode="r") as inputdata:

        lines_per_file = get_file_linecount(infile) / len(outfiles)

        file_number = 0
        line_number = 0

        outputfilehandle = get_file_handle(outfiles[file_number], mode='w')
        if headerfile:
            outputfilehandle.writelines(header)

        while True:
            line = inputdata.readline()
            line_number += 1
            if not line:
                break

            outputfilehandle.write(line)

            if line_number > ((file_number + 1) * lines_per_file):
                # if next line has the same readname as last line in file
                # then write to the same file else write to next
                line2 = inputdata.readline()
                line_number += 1
                readname2 = line2.split()[0]
                readname = line.split()[0]

                if readname2 == readname:
                    outputfilehandle.write(line2)

                outputfilehandle.close()
                outputfilehandle = open(outfiles[file_number], 'w')
                if headerfile:
                    outputfilehandle.writelines(header)

                if not readname2 == readname:
                    outputfilehandle.write(line2)

                file_number += 1

        outputfilehandle.close()
