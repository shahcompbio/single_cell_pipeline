import remixt
import remixt.seqdataio

from single_cell.utils import helpers


def seqdata_worker(chrom_seqdata, bam_file, snp_positions, chromosome,
                    bam_max_fragment_length, bam_max_soft_clipped,
                    bam_check_proper_pair):
    remixt.seqdataio.create_chromosome_seqdata(
        chrom_seqdata, bam_file, snp_positions, chromosome,
        bam_max_fragment_length, bam_max_soft_clipped,
        bam_check_proper_pair)


def create_chromosome_seqdata(seqdata, bam_file, bai_file, snp_positions, chromosomes,
                              bam_max_fragment_length, bam_max_soft_clipped,
                              bam_check_proper_pair, multiprocess=False, ncores=1):



    args = []

    for chrom in chromosomes:
        chrom_seqdata = seqdata[chrom]
        arg = (chrom_seqdata, bam_file, snp_positions, chrom,
                bam_max_fragment_length, bam_max_soft_clipped,
                bam_check_proper_pair)
        args.append(arg)

    if not multiprocess:
        for argset in args:
            seqdata_worker(*argset)
    else:
        helpers.run_in_parallel(seqdata_worker, args, ncores=ncores)

