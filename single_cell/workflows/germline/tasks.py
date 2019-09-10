import os
import shutil

import biowrappers.components.io.vcf.tasks
import pandas as pd
import vcf

import pypeliner
from pandas.api.types import CategoricalDtype

NUCLEOTIDES = ('A', 'C', 'G', 'T')


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


def run_samtools_variant_calling(
        bam_file,
        ref_genome_fasta_file,
        out_file,
        max_depth=int(1e7),
        min_bqual=0,
        min_depth=0,
        min_mqual=0,
        region=None,
        samtools_docker=None,
        vcftools_docker=None):
    mpileup_cmd = [
        'samtools',
        'mpileup',
        '-ugf', ref_genome_fasta_file,
        '-Q', min_bqual,
        '-q', min_mqual,
        bam_file
    ]

    if region is not None:
        region = '{}:{}-{}'.format(*region.split('-'))
        mpileup_cmd.extend(['-r', region])

    bcf_cmd = [
        'bcftools',
        'call',
        '-vmO', 'z',
        '-o', out_file,
    ]

    cmd = []

    cmd.extend(mpileup_cmd)
    cmd.append('|')
    cmd.extend(bcf_cmd)

    pypeliner.commandline.execute(*cmd, **samtools_docker)

    biowrappers.components.io.vcf.tasks.index_bcf(out_file, docker_config=vcftools_docker)


def annotate_normal_genotype(vcf_filename, results_filename, chromosomes):
    """ Extract het call for SNPs.
    """
    vcf_reader = vcf.Reader(filename=vcf_filename)

    hdf_store = pd.HDFStore(results_filename, 'w', complevel=9, complib='blosc')
    table_name = '/genotype'

    nucleotides = ['A', 'C', 'T', 'G']

    chunk_size = 100000

    num_rows = 0
    genotype_table = []

    def store_table(genotype_table, num_rows):
        if len(genotype_table) > 0:
            df = pd.DataFrame(
                genotype_table,
                index=range(num_rows, num_rows + len(genotype_table)),
                columns=['chrom', 'coord', 'ref', 'alt', 'is_het'])

            df['chrom'] = df['chrom'].astype(CategoricalDtype(categories=chromosomes))
            df['ref'] = df['ref'].astype(CategoricalDtype(categories=nucleotides))
            df['alt'] = df['alt'].astype(CategoricalDtype(categories=nucleotides))

            hdf_store.append(table_name, df)

            num_rows += len(genotype_table)

        return num_rows

    for record in vcf_reader:
        chrom, coord, ref, alt = record.CHROM, record.POS, record.REF, str(record.ALT[0])
        is_het = (record.samples[0].gt_alleles[0] != record.samples[0].gt_alleles[1]) * 1

        if chrom not in chromosomes:
            continue
        if ref not in nucleotides:
            continue
        if alt not in nucleotides:
            continue

        genotype_table.append([chrom, coord, ref, alt, is_het])

        if len(genotype_table) >= chunk_size:
            num_rows = store_table(genotype_table, num_rows)
            genotype_table = []

    num_rows = store_table(genotype_table, num_rows)

    hdf_store.close()
