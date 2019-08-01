'''
Created on Feb 27, 2018

@author: dgrewal
'''
import logging
import os


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
