'''
Created on Feb 27, 2018

@author: dgrewal
'''
import os
import logging
import pypeliner.commandline as cli
from single_cell.utils import helpers

def _get_header(infile):
    '''
    Extract header from the VCF file

    :param infile: input VCF file
    :return: header
    '''

    header = []
    for line in infile:
        if line.startswith('##'):
            header.append(line)
        elif line.startswith('#'):
            header.append(line)
            return header
        else:
            raise Exception('invalid header: missing #CHROM line')

    logging.getLogger("single_cell.helpers.vcfutils").warn(
        "One of the input files is empty"
    )
    return []

def concatenate_vcf(in_files, out_file, tempdir, allow_overlap=False, bcf_output=False, index_file=None, docker_image=None):
    """ Concatenate VCF files.
    :param in_files: dict with values being files to be concatenated. Files will be concatenated based on sorted order
        of keys.
    :param out_file: path where output file will be written in VCF format.
    :param allow_overlap: bool indicating whether there maybe overlapping regions in the files.
    :param index_file: path where tabix file will be written.
    """
    helpers.makedirs(tempdir)
    merged_output = os.path.join(tempdir, 'merged_vcf_files.vcf')
    sorted_output = os.path.join(tempdir, 'sorted_merged_vcf_files.vcf')

    if bcf_output:
        output_format_str = 'b'

    else:
        output_format_str = 'z'

    if isinstance(in_files, dict):
        in_files = in_files.values()

    if allow_overlap:
        cmd = ['bcftools', 'concat', '-a', '-o', merged_output]
    else:
        cmd = ['bcftools', 'concat', '-o', merged_output]

    cmd.extend(in_files)

    cli.execute(*cmd, docker_image=docker_image)

    sort_vcf(merged_output, sorted_output, docker_image=docker_image)
    bgzip_vcf(sorted_output, out_file, docker_image=docker_image)
    index_vcf(out_file, docker_image=docker_image)
    index_bcf(out_file, docker_image=docker_image)


def bgzip_vcf(infile, outfile, docker_image=None):
    cli.execute('bgzip', infile, '-c', '>', outfile, docker_image=docker_image)


def sort_vcf(in_file, out_file, docker_image=None):

    cli.execute(
        'vcf-sort',
        in_file,
        '>',
        out_file,
        docker_image=docker_image
    )


def index_vcf(vcf_file, docker_image=None):
    """ Create a tabix index for a VCF file

    :param vcf_file: Path of VCF to create index for. Should compressed by bgzip.
    :param index_file: Path of index file.

    This is meant to be used from pypeliner so it does some name mangling to add .tmp to the index file.

    """

    cli.execute(
        'tabix', '-f', '-p', 'vcf', vcf_file,
        docker_image=docker_image)

def index_bcf(in_file, docker_image):
    """ Index a VCF or BCF file with bcftools.

    :param in_file: Path of file to index.
    :param index_file: Path of index file.

    """

    cli.execute(
        'bcftools', 'index', in_file,
        docker_image=docker_image)


def filter_vcf(raw_vcf, filtered_vcf, docker_image=None):
    cmd = [
        'bcftools',
        'view',
        '-O', 'z',
        '-f', '.,PASS',
        '-o', filtered_vcf,
        raw_vcf,
    ]

    cli.execute(*cmd, docker_image=docker_image)

    index_vcf(filtered_vcf, docker_image=docker_image)
    index_bcf(filtered_vcf, docker_image=docker_image)
