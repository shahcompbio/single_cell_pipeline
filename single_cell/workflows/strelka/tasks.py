'''
Created on Nov 1, 2015

@author: Andrew Roth
'''
from __future__ import division

import csv
import math
import re
from collections import OrderedDict

from configparser import ConfigParser
import pandas as pd
import vcf

import pypeliner

FILTER_ID_BASE = 'BCNoise'
FILTER_ID_DEPTH = 'DP'
FILTER_ID_INDEL_HPOL = 'iHpol'
FILTER_ID_QSI = 'QSI_ref'
FILTER_ID_QSS = 'QSS_ref'
FILTER_ID_REPEAT = 'Repeat'
FILTER_ID_SPANNING_DELETION = 'SpanDel'


def _get_files_for_chrom(infiles, intervals, chrom):
    if not isinstance(infiles, dict):
        infiles = {ival: infiles[ival] for ival in intervals}

    outfiles = {}

    for interval in intervals:
        ival_chrom = interval.split("-")[0]

        if ival_chrom == chrom:
            outfiles[interval] = infiles[interval]

    return outfiles


def get_known_chromosome_sizes(size_file, chromosomes):
    sizes = {}

    with open(size_file, 'r') as fh:
        reader = csv.DictReader(fh, ['path', 'chrom', 'known_size', 'size'], delimiter='\t')

        for row in reader:
            if row['chrom'] not in chromosomes:
                continue

            sizes[row['chrom']] = int(row['known_size'])

    return sizes


def count_fasta_bases(ref_genome_fasta_file, out_file, docker_config):
    cmd = [
        'countFastaBases',
        ref_genome_fasta_file,
        '>',
        out_file
    ]

    pypeliner.commandline.execute(*cmd, **docker_config)


def call_somatic_variants(
        normal_bam_file,
        tumour_bam_file,
        known_sizes,
        ref_genome,
        indel_file,
        indel_window_file,
        snv_file,
        stats_file,
        region,
        docker_config,
        ncores=None,
        max_input_depth=10000,
        min_tier_one_mapq=20,
        min_tier_two_mapq=5,
        sindel_noise=0.000001,
        sindel_prior=0.000001,
        ssnv_noise=0.0000005,
        ssnv_noise_strand_bias_frac=0.5,
        ssnv_prior=0.000001):
    chrom, beg, end = re.split("-", region)

    genome_size = sum(known_sizes.values())

    beg = int(beg)
    beg = beg + 1 if beg == 0 else beg
    end = int(end)

    cmd = [
        'strelka2',
        '-bam-file', normal_bam_file,
        '--tumor-bam-file', tumour_bam_file,

        '-samtools-reference', ref_genome,

        '--somatic-indel-file', indel_file,
        '--somatic-snv-file', snv_file,
        '--variant-window-flank-file', 50, indel_window_file,

        '-bam-seq-name', chrom,
        '-report-range-begin', beg,
        '-report-range-end', end,

        '-clobber',
        '-filter-unanchored',
        '-genome-size', genome_size,
        '-indel-nonsite-match-prob', 0.5,
        '-max-indel-size', 50,
        '-max-window-mismatch', 3, 20,
        '-min-paired-align-score', min_tier_one_mapq,
        '-min-single-align-score', 10,
        '-min-qscore', 0,
        '-print-used-allele-counts',

        '--max-input-depth', max_input_depth,
        '--min-contig-open-end-support', 35,
        '--report-file', stats_file,
        '--shared-site-error-rate', ssnv_noise,

        '--shared-site-error-strand-bias-fraction', ssnv_noise_strand_bias_frac,
        '--somatic-indel-rate', sindel_prior,
        '--shared-indel-error-rate', sindel_noise,
        '--somatic-snv-rate', ssnv_prior,
        '--tier2-include-anomalous',
        '--tier2-include-singleton',
        '--tier2-indel-nonsite-match-prob', 0.25,
        '--tier2-min-paired-align-score', min_tier_two_mapq,
        '--tier2-min-single-align-score', min_tier_two_mapq,
        '--tier2-mismatch-density-filter-count', 10,
        '--tier2-no-filter-unanchored',
        '--tier2-single-align-score-rescue-mode'
    ]

    pypeliner.commandline.execute(*cmd, **docker_config)


# =======================================================================================================================
# SNV filtering
# =======================================================================================================================


def filter_snv_file_list(
        in_files,
        stats_files,
        out_file,
        chrom,
        known_chrom_size,
        intervals,
        depth_filter_multiple=3.0,
        max_filtered_basecall_frac=0.4,
        max_spanning_deletion_frac=0.75,
        quality_lower_bound=15,
        use_depth_filter=True):
    known_chrom_size = known_chrom_size[chrom]

    in_files = _get_files_for_chrom(in_files, intervals, chrom)
    stats_files = _get_files_for_chrom(stats_files, intervals, chrom)

    max_normal_coverage = _get_max_normal_coverage(chrom, depth_filter_multiple, known_chrom_size, stats_files)

    writer = None

    with open(out_file, 'wt') as out_fh:
        for key in sorted(in_files):
            reader = vcf.Reader(filename=in_files[key])

            if writer is None:
                # Add filters to header
                if use_depth_filter:
                    reader.filters[FILTER_ID_DEPTH] = vcf.parser._Filter(
                        id=FILTER_ID_DEPTH,
                        desc='Greater than {0}x chromosomal mean depth in Normal sample'.format(depth_filter_multiple)
                    )

                reader.filters[FILTER_ID_BASE] = vcf.parser._Filter(
                    id=FILTER_ID_BASE,
                    desc='Fraction of basecalls filtered at this site in either sample is at or above {0}'.format(
                        max_filtered_basecall_frac)
                )

                reader.filters[FILTER_ID_SPANNING_DELETION] = vcf.parser._Filter(
                    id=FILTER_ID_SPANNING_DELETION,
                    desc='Fraction of reads crossing site with spanning deletions in either sample exceeeds {0}'.format(
                        max_spanning_deletion_frac)
                )

                reader.filters[FILTER_ID_QSS] = vcf.parser._Filter(
                    id=FILTER_ID_QSS,
                    desc='Normal sample is not homozygous ref or ssnv Q-score < {0}, ie calls with NT!=ref or QSS_NT < {0}'.format(
                        quality_lower_bound)
                )

                writer = vcf.Writer(out_fh, reader)

            for record in reader:
                normal = record.genotype('NORMAL')

                tumour = record.genotype('TUMOR')

                # Normal depth filter
                if use_depth_filter and (normal.data.DP > max_normal_coverage):
                    record.add_filter(FILTER_ID_DEPTH)

                # Filtered basecall fraction
                normal_filtered_base_call_fraction = _get_filtered_base_call_fraction(normal.data)

                tumour_filtered_base_call_fraction = _get_filtered_base_call_fraction(tumour.data)

                if (normal_filtered_base_call_fraction >= max_filtered_basecall_frac) or \
                        (tumour_filtered_base_call_fraction >= max_filtered_basecall_frac):
                    record.add_filter(FILTER_ID_BASE)

                # Spanning deletion fraction
                normal_spanning_deletion_fraction = _get_spanning_deletion_fraction(normal.data)

                tumour_spanning_deletion_fraction = _get_spanning_deletion_fraction(tumour.data)

                if (normal_spanning_deletion_fraction > max_spanning_deletion_frac) or \
                        (tumour_spanning_deletion_fraction > max_spanning_deletion_frac):
                    record.add_filter(FILTER_ID_SPANNING_DELETION)

                # Q-val filter
                if (record.INFO['NT'] != 'ref') or (record.INFO['QSS_NT'] < quality_lower_bound):
                    record.add_filter(FILTER_ID_QSS)

                writer.write_record(record)

        writer.close()


def _get_max_normal_coverage(chrom, depth_filter_multiple, known_chrom_size, stats_files):
    normal_coverage = _get_normal_coverage(stats_files.values())

    normal_mean_coverage = normal_coverage / known_chrom_size

    max_normal_coverage = normal_mean_coverage * depth_filter_multiple

    return max_normal_coverage


def _get_normal_coverage(stats_files):
    total_coverage = 0

    for file_name in stats_files:
        mean_matcher = re.compile('mean:\s(.*?)\s')

        sample_size_matcher = re.compile('sample_size:\s(.*?)\s')

        with open(file_name) as fh:
            for line in fh:
                if line.startswith('NORMAL_NO_REF_N_COVERAGE '):
                    mean = float(mean_matcher.search(line).group(1))

                    sample_size = float(sample_size_matcher.search(line).group(1))

                    if math.isnan(mean) or math.isnan(sample_size):
                        continue

                    total_coverage += mean * sample_size

    return total_coverage


def _get_filtered_base_call_fraction(data):
    frac = 0

    if data.DP > 0:
        frac = data.FDP / data.DP

    return frac


def _get_spanning_deletion_fraction(data):
    total = data.DP + data.SDP

    frac = 0

    if total > 0:
        frac = data.SDP / total

    return frac


# =======================================================================================================================
# Indel filtering
# =======================================================================================================================


def filter_indel_file_list(
        vcf_files,
        stats_files,
        window_files,
        out_file,
        chrom,
        known_chrom_size,
        intervals,
        depth_filter_multiple=3.0,
        max_int_hpol_length=14,
        max_ref_repeat=8,
        max_window_filtered_basecall_frac=0.3,
        quality_lower_bound=30,
        use_depth_filter=True):
    window_cols = (
        'chrom',
        'coord',
        'normal_window_used',
        'normal_window_filtered',
        'normal_window_submap',
        'tumour_window_used',
        'tumour_window_filtered',
        'tumour_window_submap'
    )

    known_chrom_size = known_chrom_size[chrom]

    vcf_files = _get_files_for_chrom(vcf_files, intervals, chrom)
    stats_files = _get_files_for_chrom(stats_files, intervals, chrom)
    window_files = _get_files_for_chrom(window_files, intervals, chrom)

    max_normal_coverage = _get_max_normal_coverage(chrom, depth_filter_multiple, known_chrom_size, stats_files)

    writer = None

    with open(out_file, 'wt') as out_fh:
        for key in sorted(vcf_files):
            window = pd.read_csv(
                window_files[key],
                comment='#',
                converters={'chrom': str},
                header=None,
                names=window_cols,
                sep='\t')

            reader = vcf.Reader(filename=vcf_files[key])

            if writer is None:
                # Add format to header
                reader.formats['DP50'] = vcf.parser._Format(
                    id='DP50',
                    num=1,
                    type='Float',
                    desc='Average tier1 read depth within 50 bases'
                )

                reader.formats['FDP50'] = vcf.parser._Format(
                    id='FDP50',
                    num=1,
                    type='Float',
                    desc='Average tier1 number of basecalls filtered from original read depth within 50 bases'
                )

                reader.formats['SUBDP50'] = vcf.parser._Format(
                    id='SUBDP50',
                    num=1,
                    type='Float',
                    desc='Average number of reads below tier1 mapping quality threshold aligned across sites within 50 bases'
                )

                # Add filters to header
                if use_depth_filter:
                    reader.filters[FILTER_ID_DEPTH] = vcf.parser._Filter(
                        id=FILTER_ID_DEPTH,
                        desc='Greater than {0}x chromosomal mean depth in Normal sample'.format(depth_filter_multiple)
                    )

                reader.filters[FILTER_ID_REPEAT] = vcf.parser._Filter(
                    id=FILTER_ID_REPEAT,
                    desc='Sequence repeat of more than {0}x in the reference sequence'.format(max_ref_repeat)
                )

                reader.filters[FILTER_ID_INDEL_HPOL] = vcf.parser._Filter(
                    id=FILTER_ID_INDEL_HPOL,
                    desc='Indel overlaps an interrupted homopolymer longer than {0}x in the reference sequence'.format(
                        max_int_hpol_length)
                )

                reader.filters[FILTER_ID_BASE] = vcf.parser._Filter(
                    id=FILTER_ID_BASE,
                    desc='Average fraction of filtered basecalls within 50 bases of the indel exceeds {0}'.format(
                        max_window_filtered_basecall_frac)
                )

                reader.filters[FILTER_ID_QSI] = vcf.parser._Filter(
                    id=FILTER_ID_QSI,
                    desc='Normal sample is not homozygous ref or sindel Q-score < {0}, ie calls with NT!=ref or QSI_NT < {0}'.format(
                        quality_lower_bound)
                )

                writer = vcf.Writer(out_fh, reader)

            for record in reader:
                window_row = window.loc[
                    (window['chrom'] == str(record.CHROM)) & (window['coord'] == record.POS)].iloc[0]

                normal = record.genotype('NORMAL')

                tumour = record.genotype('TUMOR')

                normal_data = normal.data._asdict()

                tumour_data = tumour.data._asdict()

                # Add window data to vcf record
                record.add_format('DP50')

                record.add_format('FDP50')

                record.add_format('SUBDP50')

                normal_data['DP50'] = window_row['normal_window_used'] + window_row['normal_window_filtered']

                normal_data['FDP50'] = window_row['normal_window_filtered']

                normal_data['SUBDP50'] = window_row['normal_window_submap']

                tumour_data['DP50'] = window_row['tumour_window_used'] + window_row['tumour_window_filtered']

                tumour_data['FDP50'] = window_row['tumour_window_filtered']

                tumour_data['SUBDP50'] = window_row['tumour_window_submap']

                normal.data = _convert_dict_to_call(normal_data)

                tumour.data = _convert_dict_to_call(tumour_data)

                # Add filters

                # Normal depth filter
                if use_depth_filter and (normal.data.DP > max_normal_coverage):
                    record.add_filter(FILTER_ID_DEPTH)

                # Ref repeat
                if 'RC' in record.INFO:
                    ref_repeat = record.INFO['RC']

                    if ref_repeat > max_ref_repeat:
                        record.add_filter(FILTER_ID_REPEAT)

                # Indel homopolymer
                if 'IHP' in record.INFO:
                    indel_homopolymer = record.INFO['IHP']

                    if indel_homopolymer > max_int_hpol_length:
                        record.add_filter(FILTER_ID_INDEL_HPOL)

                # Base filter
                normal_filtered_base_call_fraction = _get_filtered_base_call_fraction_indel(normal.data)

                tumour_filtered_base_call_fraction = _get_filtered_base_call_fraction_indel(tumour.data)

                if (normal_filtered_base_call_fraction >= max_window_filtered_basecall_frac) or \
                        (tumour_filtered_base_call_fraction >= max_window_filtered_basecall_frac):
                    record.add_filter(FILTER_ID_BASE)

                # Q-val filter
                if (record.INFO['NT'] != 'ref') or (record.INFO['QSI_NT'] < quality_lower_bound):
                    record.add_filter(FILTER_ID_QSI)

                writer.write_record(record)

        if writer:
            writer.close()


def _convert_dict_to_call(data_dict):
    call_data_class = vcf.model.make_calldata_tuple(data_dict.keys())

    return call_data_class(**data_dict)


def _get_filtered_base_call_fraction_indel(data):
    frac = 0

    if data.DP50 > 0:
        frac = data.FDP50 / data.DP50

    return frac


# =======================================================================================================================
# Write config file for make style strelka
# =======================================================================================================================


def configure_run(normal_bam_file,
                  tumour_bam_file,
                  ref_genome_fasta_file,
                  ref_genome_base_counts_file,
                  out_file,
                  bin_size=int(1e7),
                  depth_filter_multiple=3.0,
                  extra_strelka_arguments='',
                  indel_max_int_hpol_length=14,
                  indel_max_ref_repeat=8,
                  indel_max_window_filtered_basecall_frac=0.3,
                  max_input_depth=10000,
                  min_tier_one_mapq=20,
                  min_tier_two_mapq=5,
                  sindel_noise=0.000001,
                  sindel_prior=0.000001,
                  sindel_quality_lower_bound=30,
                  snv_max_filtered_basecall_frac=0.4,
                  snv_max_spanning_deletion_frac=0.75,
                  ssnv_noise=0.0000005,
                  ssnv_noise_strand_bias_frac=0.5,
                  ssnv_prior=0.000001,
                  ssnv_quality_lower_bound=15,
                  skip_depth_filters=False,
                  write_realigned_bam=False):
    parser = ConfigParser.ConfigParser()

    parser.add_section('user')

    parser.set('user', 'binSize', bin_size)
    parser.set('user', 'depthFilterMultiple', depth_filter_multiple)
    parser.set('user', 'extraStrelkaArguments', extra_strelka_arguments)
    parser.set('user', 'indelMaxIntHpolLength', indel_max_int_hpol_length)
    parser.set('user', 'indelMaxRefRepeat', indel_max_ref_repeat)
    parser.set('user', 'indelMaxWindowFilteredBasecallFrac', indel_max_window_filtered_basecall_frac)
    parser.set('user', 'isSkipDepthFilters', int(skip_depth_filters))
    parser.set('user', 'isWriteRealignedBam', int(write_realigned_bam))
    parser.set('user', 'maxInputDepth', max_input_depth)
    parser.set('user', 'minTier1Mapq', min_tier_one_mapq)
    parser.set('user', 'minTier2Mapq', min_tier_two_mapq)
    parser.set('user', 'sindelNoise', sindel_noise)
    parser.set('user', 'sindelPrior', sindel_prior)
    parser.set('user', 'sindelQuality_LowerBound', sindel_quality_lower_bound)
    parser.set('user', 'snvMaxFilteredBasecallFrac', snv_max_filtered_basecall_frac)
    parser.set('user', 'snvMaxSpanningDeletionFrac', snv_max_spanning_deletion_frac)
    parser.set('user', 'ssnvNoise', ssnv_noise)
    parser.set('user', 'ssnvNoiseStrandBiasFrac', ssnv_noise_strand_bias_frac)
    parser.set('user', 'ssnvPrior', ssnv_prior)
    parser.set('user', 'ssnvQuality_LowerBound', ssnv_quality_lower_bound)

    parser.add_section('derived')

    chrom_order = []

    known_genome_size = 0

    with open(ref_genome_base_counts_file) as fh:
        reader = csv.DictReader(fh, ['path', 'chrom', 'known_size', 'size'], delimiter='\t')

        for row in reader:
            chrom_order.append(row['chrom'])

            parser.set('derived', 'chrom_{0}_knownSize'.format(row['chrom']), row['known_size'])

            parser.set('derived', 'chrom_{0}_size'.format(row['chrom']), row['size'])

            known_genome_size += int(row['known_size'])

    parser.set('derived', 'chromOrder', '\t'.join(chrom_order))

    parser.set('derived', 'knownGenomeSize', known_genome_size)

    parser.set('derived', 'normalBam', normal_bam_file)

    parser.set('derived', 'tumorBam', tumour_bam_file)

    parser.set('derived', 'refFile', ref_genome_fasta_file)

    with open(out_file, 'wt') as out_fh:
        parser.write(out_fh)


# =======================================================================================================================
# Dump to table
# =======================================================================================================================


def convert_vcf_to_hdf5(in_file, out_file, data_type='snv', table_name=None):
    out_store = pd.HDFStore(out_file, 'w', complevel=9, complib='blosc')

    reader = vcf.Reader(filename=in_file)

    if table_name is None:
        table_name = 'strelka_{0}'.format(data_type)

    for record in reader:
        row = OrderedDict()

        row['chrom'] = str(record.CHROM)

        row['coord'] = int(record.POS)

        row['ref_base'] = str(record.REF)

        row['alt_base'] = str(record.ALT[0])

        if data_type == 'snv':
            row['qual'] = float(record.INFO['QSS'])

        elif data_type == 'indel':
            row['qual'] = float(record.INFO['QSI'])

        else:
            raise Exception('Unknown data type {0}'.format(data_type))

        row = pd.DataFrame([pd.Series(row)])

        out_store.append(table_name, row)

    out_store.close()
