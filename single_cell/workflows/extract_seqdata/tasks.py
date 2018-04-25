import remixt
import remixt.seqdataio


def create_chromosome_seqdata(seqdata, bam_file, bai_file, snp_positions, chromosomes,
                              bam_max_fragment_length, bam_max_soft_clipped,
                              bam_check_proper_pair):

    for chromosome in chromosomes:

        chrom_seqdata = seqdata[chromosome]

        remixt.seqdataio.create_chromosome_seqdata(
            chrom_seqdata, bam_file, snp_positions, chromosome,
            bam_max_fragment_length, bam_max_soft_clipped,
            bam_check_proper_pair)
