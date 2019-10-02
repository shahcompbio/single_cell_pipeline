'''
Created on Feb 27, 2018

@author: dgrewal
'''
import logging
import os

import biowrappers.components.io.vcf.tasks as vcf_tasks
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


def concatenate_vcf(infiles, outfile):
    '''
    Concatenate VCF files

    :param infiles: dictionary of input VCF files to be concatenated
    :param outfile: output VCF file
    '''

    with open(outfile, 'w') as ofile:
        header = None

        for _, ifile in infiles.items():

            if os.path.getsize(ifile) == 0:
                logging.getLogger("single_cell.helpers.vcfutils").warn(
                    'input file {} is empty'.format(ifile)
                )
                continue

            with open(ifile) as f:

                if not header:
                    header = _get_header(f)

                    for line in header:
                        ofile.write(line)
                else:
                    if not _get_header(f) == header:
                        logging.getLogger("single_cell.helpers.vcfutils").warn(
                            'merging vcf files with mismatching headers'
                        )

                for l in f:
                    ofile.write(l)


def merge_vcf(infiles, outfile, tempdir, docker_image=None):
    vcf_files = []
    for infile in infiles:
        if isinstance(infile, str):
            vcf_files.append(infile)
        elif isinstance(infile, dict):
            vcf_files.extend(list(infile.values()))
        elif isinstance(infile, (list, tuple)):
            vcf_files.extend(list(infile))
        else:
            raise Exception("unknown data type")

    helpers.makedirs(tempdir)
    temp_output = os.path.join(tempdir, 'merged.vcf')

    vcf_tasks.merge_vcfs(vcf_files, temp_output)

    vcf_tasks.finalise_vcf(temp_output, outfile, docker_config={'docker_image': docker_image})
