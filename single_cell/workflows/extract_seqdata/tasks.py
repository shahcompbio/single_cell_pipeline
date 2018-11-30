import remixt
import remixt.seqdataio
import remixt.config

from single_cell.utils import helpers


def seqdata_worker(chrom_seqdata, bam_file, snp_positions, chromosome,
                    bam_max_fragment_length, bam_max_soft_clipped,
                    bam_check_proper_pair):

    remixt.seqdataio.create_chromosome_seqdata(
        chrom_seqdata, bam_file, snp_positions, chromosome,
        bam_max_fragment_length, bam_max_soft_clipped,
        bam_check_proper_pair)


def create_chromosome_seqdata(seqdata, bam_file, config, ref_data_dir,
                              multiprocess=False, ncores=1, chromosomes=None):

    if not chromosomes:
        chromosomes = remixt.config.get_chromosomes(config, ref_data_dir)

    snp_positions_filename = remixt.config.get_filename(config, ref_data_dir, 'snp_positions')

    bam_max_fragment_length = remixt.config.get_param(config, 'bam_max_fragment_length')
    bam_max_soft_clipped = remixt.config.get_param(config, 'bam_max_soft_clipped')
    bam_check_proper_pair = remixt.config.get_param(config, 'bam_check_proper_pair')

    args = []

    for chrom in chromosomes:
        chrom_seqdata = seqdata[chrom]
        arg = (
            chrom_seqdata, bam_file, snp_positions_filename, chrom,
            bam_max_fragment_length, bam_max_soft_clipped,
            bam_check_proper_pair)
        args.append(arg)

    if not multiprocess:
        for argset in args:
            seqdata_worker(*argset)
    else:
        helpers.run_in_parallel(seqdata_worker, args, ncores=ncores)
