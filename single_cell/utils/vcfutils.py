'''
Created on Feb 27, 2018

@author: dgrewal
'''
import itertools
import logging
import os

import biowrappers.components.io.vcf.tasks as vcf_tasks
import vcf
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


def merge_vcf(infiles, outfile, tempdir):
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

    vcf_tasks.finalise_vcf(temp_output, outfile)


def split_vcf(in_file, out_files, lines_per_file):
    """ Split a VCF file into smaller files.

    :param in_file: Path of VCF file to split.

    :param out_files: Callback function which supplies file name given index of split.

    :param lines_per_file: Maximum number of lines to be written per file.

     """

    def line_group(_, line_idx=itertools.count()):
        return int(next(line_idx) / lines_per_file)

    reader = vcf.Reader(filename=in_file)

    for file_idx, records in itertools.groupby(reader, key=line_group):
        file_name = out_files[file_idx]

        with open(file_name, 'w') as out_fh:
            writer = vcf.Writer(out_fh, reader)

            for record in records:
                writer.write_record(record)

            writer.close()
