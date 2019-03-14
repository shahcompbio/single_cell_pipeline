import os
import remixt
import remixt.seqdataio
import remixt.config

from single_cell.utils import helpers


def create_chromosome_seqdata(seqdata, bam_file, tempdir, config, ref_data_dir,
                              chromosomes=None):
    helpers.makedirs(tempdir)
    if not chromosomes:
        chromosomes = remixt.config.get_chromosomes(config, ref_data_dir)

    snp_positions_filename = remixt.config.get_filename(config, ref_data_dir, 'snp_positions')

    all_seqdata = {}

    bam_max_fragment_length = remixt.config.get_param(config, 'bam_max_fragment_length')
    bam_max_soft_clipped = remixt.config.get_param(config, 'bam_max_soft_clipped')
    bam_check_proper_pair = remixt.config.get_param(config, 'bam_check_proper_pair')

    for chrom in chromosomes:
        chrom_seqdata = os.path.join(tempdir, "{}_seqdata.h5".format(chrom))
        all_seqdata[chrom] = chrom_seqdata

        remixt.seqdataio.create_chromosome_seqdata(
            chrom_seqdata, bam_file, snp_positions_filename,
            chrom, bam_max_fragment_length, bam_max_soft_clipped,
            bam_check_proper_pair)

    remixt.seqdataio.merge_seqdata(seqdata, all_seqdata)