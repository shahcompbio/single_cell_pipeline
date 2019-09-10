import logging

import numpy as np
import pandas as pd
import remixt
from .scripts import merge_wigs
from single_cell.utils import csvutils
from pandas.api.types import CategoricalDtype


def merge_overlapping_seqdata(outfile, infiles, chromosomes):
    out_store = pd.HDFStore(outfile, 'w', complevel=9, complib='blosc')

    index_offsets = pd.Series(0, index=chromosomes, dtype=np.int64)

    for _id, infile in infiles.items():
        store = pd.HDFStore(infile)
        tables = store.keys()

        for chromosome in chromosomes:
            allele_table = '/alleles/chromosome_{}'.format(chromosome)
            fragment_table = '/fragments/chromosome_{}'.format(chromosome)

            if allele_table not in tables:
                logging.getLogger("single_cell.titan").warn(
                    "missing table {}".format(allele_table)
                )
                continue

            if fragment_table not in tables:
                logging.getLogger("single_cell.titan").warn(
                    "missing table {}".format(fragment_table)
                )
                continue

            alleles = store[allele_table]
            fragments = store[fragment_table]

            alleles['fragment_id'] = alleles['fragment_id'].astype(np.int64)
            fragments['fragment_id'] = fragments['fragment_id'].astype(np.int64)

            alleles['fragment_id'] += index_offsets[chromosome]
            fragments['fragment_id'] += index_offsets[chromosome]

            index_offsets[chromosome] = max(alleles['fragment_id'].max(), fragments['fragment_id'].max()) + 1

            out_store.append('/alleles/chromosome_{}'.format(chromosome), alleles)
            out_store.append('/fragments/chromosome_{}'.format(chromosome), fragments)

        store.close()

    out_store.close()


def create_chromosome_seqdata(seqdata, bam_file, snp_positions, chromosomes,
                              bam_max_fragment_length, bam_max_soft_clipped,
                              bam_check_proper_pair):
    for chromosome in chromosomes:
        chrom_seqdata = seqdata[chromosome]

        remixt.seqdataio.create_chromosome_seqdata(
            chrom_seqdata, bam_file, snp_positions, chromosome,
            bam_max_fragment_length, bam_max_soft_clipped,
            bam_check_proper_pair)


def merge_wig_files(input_wigs, output):
    input_wigs = input_wigs.values()
    merge_wigs.main(input_wigs, output)


def merge_het_positions(input_csvs, output):
    csvutils.merge_csv(
        input_csvs, output, how="outer", on=[
            "chromosome", "position"], sep="\t")


def merge_tumour_alleles(input_csvs, output):
    data = {}

    for _, infile in input_csvs.items():

        with open(infile) as reader:
            for line in reader:
                line = line.strip().split()
                chrom, pos, _, ref, _, alt = line
                ref = int(ref)
                alt = int(alt)

                if (chrom, pos) not in data:
                    data[(chrom, pos)] = (ref, alt)
                else:
                    oldref, oldalt = data[(chrom, pos)]
                    data[(chrom, pos)] = (oldref + ref, oldalt + alt)

    with open(output, "w") as writer:
        for (chrom, pos), (ref, alt) in data.items():
            writer.write("{}\t{}\tA\t{}\tT\t{}\n".format(chrom, pos, ref, alt))


def concat_tumour_alleles(input_csvs, output_filename, chromosomes):
    store = pd.HDFStore(output_filename, 'w', complevel=9, complib='blosc')

    for cell_id, input_filename in input_csvs.items():
        cell_data = pd.read_csv(input_filename, sep='\t', header=None,
                                names=['chromosome', 'coord', 'ref', 'ref_counts', 'alt', 'alt_counts'])
        cell_data = cell_data.drop(['ref', 'alt'], axis=1)
        cell_data['chromosome'] = cell_data['chromosome'].astype(CategoricalDtype(categories=chromosomes))
        cell_data['cell_id'] = cell_id
        cell_data['cell_id'] = cell_data['cell_id'].astype(CategoricalDtype(categories=input_csvs.keys()))

        store.append('/allele_counts', cell_data, data_columns=['cell_id', 'chromosome'])

    store.close()
