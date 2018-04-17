from biowrappers.components.copy_number_calling import titan
import pandas as pd
import itertools

from scripts import merge_wigs

from single_cell.utils import csvutils
import remixt


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

    for _, infile in input_csvs.iteritems():

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
        for (chrom, pos), (ref, alt) in data.iteritems():
            writer.write("{}\t{}\tA\t{}\tT\t{}\n".format(chrom, pos, ref, alt))
