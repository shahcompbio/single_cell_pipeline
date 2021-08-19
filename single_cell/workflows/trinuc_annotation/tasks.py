import pandas as pd
import pysam
import vcf

from single_cell.utils import csvutils
from single_cell.workflows.trinuc_annotation.dtypes import dtypes


def get_tri_nucelotide_context(ref_genome_fasta_file, vcf_file, out_file):
    vcf_reader = vcf.Reader(filename=vcf_file)

    fasta_reader = pysam.Fastafile(ref_genome_fasta_file)

    data = []

    for record in vcf_reader:
        chrom = record.CHROM

        coord = record.POS

        tri_nucleotide_context = fasta_reader.fetch(chrom, coord - 2, coord + 1)

        data.append({'chrom': record.CHROM, 'coord': record.POS, 'tri_nucleotide_context': tri_nucleotide_context})

    data = pd.DataFrame(data)

    csvutils.write_dataframe_to_csv_and_yaml(data, out_file, dtypes())
