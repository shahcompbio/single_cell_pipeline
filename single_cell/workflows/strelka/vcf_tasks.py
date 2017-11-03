'''
Created on Oct 31, 2015

@author: Andrew Roth
'''
# from pypeliner.workflow import Workflow

import os
import itertools
import pandas as pd
import pypeliner
import shutil
# import time
import vcf

from components_utils import flatten_input

# from ._merge import merge_vcfs


def compress_vcf(in_file, out_file, index_file=None):
    """ Compress a VCF file using bgzip.

    :param in_file: Path of uncompressed VCF file.
    :param out_file: Path were compressed VCF file will be written.
    """
    pypeliner.commandline.execute('bgzip', '-c', in_file, '>', out_file)
    if index_file is not None:
        index_vcf(out_file, index_file=index_file)


def filter_vcf(in_file, out_file):
    """ Filter a VCF for records with no filters set.

    :param in_file: Path of VCF file to filter.

    :param out_file: Path where filtered VCF file will be written.

    Note that records with the filter `PASS` will not be removed.

    """

    reader = vcf.Reader(filename=in_file)

    with open(out_file, 'wb') as out_fh:
        writer = vcf.Writer(out_fh, reader)

        for record in reader:
            if (record.FILTER is None) or (len(record.FILTER) == 0):
                writer.write_record(record)

        writer.close()


def _rename_index(in_file, index_suffix):
    if in_file.endswith('.tmp'):
        index_file = in_file[:-4] + index_suffix

        try:
            os.remove(index_file)
        except:
            pass

        shutil.move(in_file + index_suffix, index_file)


def index_bcf(in_file, index_file=None):
    """ Index a VCF or BCF file with bcftools.

    :param in_file: Path of file to index.
    :param index_file: Path of index file.

    """

    pypeliner.commandline.execute('bcftools', 'index', in_file)

    if index_file is None:
        _rename_index(in_file, '.csi')

    else:
        shutil.move(in_file + '.csi', index_file)


def finalise_vcf(in_file, compressed_file):
    """ Compress a VCF using bgzip and create index.

    :param in_file: Path of file to compressed and index.

    :param out_file: Path where compressed file will be written. Index file will written to `out_file` + `.tbi` and `out_file` + `.csi` and .

    """

    uncompressed_file = compressed_file + '.uncompressed'
    pypeliner.commandline.execute('vcf-sort', in_file, '>', uncompressed_file)
    pypeliner.commandline.execute('bgzip', uncompressed_file, '-c', '>', compressed_file)
    os.remove(uncompressed_file)

    index_bcf(compressed_file)
    index_vcf(compressed_file)


def index_vcf(vcf_file, index_file=None):
    """ Create a tabix index for a VCF file

    :param vcf_file: Path of VCF to create index for. Should compressed by bgzip.
    :param index_file: Path of index file.

    This is meant to be used from pypeliner so it does some name mangling to add .tmp to the index file.

    """

    pypeliner.commandline.execute('tabix', '-f', '-p', 'vcf', vcf_file)

    if index_file is None:
        _rename_index(vcf_file, '.tbi')

    else:
        shutil.move(vcf_file + '.tbi', index_file)


def concatenate_vcf(in_files, out_file, allow_overlap=False, bcf_index_file=None, vcf_index_file=None):
    """ Fast concatenation of VCF file using `bcftools`.

    :param in_files: dict with values being files to be concatenated. Files will be concatenated based on sorted order of keys.

    :param out_file: path where output file will be written in VCF format.

    """
    if allow_overlap:
        cmd = ['bcftools', 'concat', '-a', '-O', 'z', '-o', out_file]
    else:
        cmd = ['bcftools', 'concat', '-O', 'z', '-o', out_file]

    cmd += flatten_input(in_files)

    pypeliner.commandline.execute(*cmd)

    index_vcf(out_file, index_file=vcf_index_file)
    index_bcf(out_file, index_file=bcf_index_file)


def concatenate_bcf(in_files, out_file):
    """ Fast concatenation of BCF file using `bcftools`.

    :param in_files: dict with values being files to be concatenated. Files will be concatenated based on sorted order of keys.

    :param out_file: path where output file will be written in VCF format.

    """

    cmd = ['bcftools', 'concat', '-a', '-O', 'b', '-o', out_file]
    cmd += flatten_input(in_files)

    pypeliner.commandline.execute(*cmd)

    index_vcf(out_file)
    index_bcf(out_file)


def extract_variant_type(in_file, out_file, variant_type):
    """ Extract a specific type of variant from a vcf file.

    :param in_file: input vcf file
    :param out_file: output filtered vcf file
    :param variant_type: one of snps or indels

    """

    pypeliner.commandline.execute('bcftools', 'view', '-v', variant_type, '-O', 'z', '-o', out_file, in_file)

    index_vcf(out_file)
    index_bcf(out_file)


def split_vcf(in_file, out_file_callback, lines_per_file):
    """ Split a VCF file into smaller files.

    :param in_file: Path of VCF file to split.

    :param out_file_callback: Callback function which supplies file name given index of split.

    :param lines_per_file: Maximum number of lines to be written per file.

     """

    def line_group(line, line_idx=itertools.count()):
        return int(next(line_idx) / lines_per_file)

    reader = vcf.Reader(filename=in_file)

    for file_idx, records in itertools.groupby(reader, key=line_group):
        file_name = out_file_callback(file_idx)

        with open(file_name, 'w') as out_fh:
            writer = vcf.Writer(out_fh, reader)

            for record in records:
                writer.write_record(record)

            writer.close()


def convert_vcf_to_hdf5(in_file, out_file, table_name, score_callback=None):

    def line_group(line, line_idx=itertools.count()):
        return int(next(line_idx) / chunk_size)

    chunk_size = 1000

    #===================================================================================================================
    # find all entries in categories
    #===================================================================================================================
    reader = vcf.Reader(filename=in_file)

    chrom_categories = set()

    ref_categories = set()

    alt_categories = set()

    for record in reader:

        chrom_categories.add(str(record.CHROM))

        ref_categories.add(str(record.REF))

        for alt in record.ALT:

            alt_categories.add(str(alt))

    chrom_categories = sorted(chrom_categories)

    ref_categories = sorted(ref_categories)

    alt_categories = sorted(alt_categories)

    #===================================================================================================================
    # convert
    #===================================================================================================================

    # reopen reader to restart iter
    reader = vcf.Reader(filename=in_file)

    hdf_store = pd.HDFStore(out_file, 'w', complevel=9, complib='blosc')

    for file_idx, records in itertools.groupby(reader, key=line_group):

        df = []

        for record in records:
            if score_callback is not None:
                score = score_callback(record)

            else:
                score = record.QUAL

            for alt in record.ALT:
                row = {
                    'chrom': record.CHROM,
                    'coord': record.POS,
                    'ref': str(record.REF),
                    'alt': str(alt),
                    'score': score
                }

                df.append(row)

        beg = file_idx * chunk_size

        end = beg + len(df)

        df = pd.DataFrame(df, index=range(beg, end))

        df['chrom'] = df['chrom'].astype('category', categories=chrom_categories)

        df['alt'] = df['alt'].astype('category', categories=alt_categories)

        df['ref'] = df['ref'].astype('category', categories=ref_categories)

        df = df[['chrom', 'coord', 'ref', 'alt', 'score']]

        hdf_store.append(table_name, df)

    hdf_store.close()


def sort_vcf(in_file, out_file):

    pypeliner.commandline.execute(
        'vcf-sort',
        in_file,
        '>',
        out_file
    )
