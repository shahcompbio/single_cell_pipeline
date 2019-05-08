'''
Created on Nov 1, 2015

@author: Andrew Roth
'''
from __future__ import division

import ConfigParser
import csv
import math
import os
import re
import shutil
from collections import OrderedDict

import numpy as np
import pandas as pd
import pypeliner
import vcf


def get_chromosome_depth(chrom, bam_file, ref_genome, out_file, docker_image=None):
    chrom = chrom.split('-')[0]

    cmd = [
        'GetChromDepth',
        '--align-file', bam_file,
        '--chrom', chrom,
        '--output-file', out_file,
        # '--ref', ref_genome,
    ]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def merge_chromosome_depths_weighted(infiles, outfile):
    data = {}

    for region, infile in infiles.iteritems():
        size = int(region.split('-')[2]) - int(region.split('-')[1])
        with open(infile) as indata:
            depthdata = indata.readline()
            chrom, depth = depthdata.strip().split()
            if chrom not in data:
                data[chrom] = []
            data[chrom].append((float(depth)*size, size))

    with open(outfile, 'w') as output:
        for chrom, depths in data.iteritems():
            total = sum([v[0] for v in depths])
            size = sum([v[1] for v in depths])
            output.write('{}\t{}\n'.format(chrom, total/size))



def merge_chromosome_depths_plain(infiles, outfile):
    data = {}

    if isinstance(infiles, dict):
        infiles = infiles.values()

    for infile in infiles:
        with open(infile) as indata:
            depthdata = indata.readline()
            chrom, depth = depthdata.strip().split()
            if chrom not in data:
                data[chrom] = []
            data[chrom].append(float(depth))

    with open(outfile, 'w') as output:
        for chrom, depths in data.iteritems():
            output.write('{}\t{}\n'.format(chrom, np.mean(depths)))


def merge_chromosome_depths(infiles, outfile):
    if isinstance(infiles, dict):
        merge_chromosome_depths_weighted(infiles, outfile)
    else:
        merge_chromosome_depths_plain(infiles, outfile)


def call_genome_segment(
        chrom_depth_file,
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        indel_file,
        snv_file,
        tmp_dir,
        region,
        known_sizes,
        is_exome=False,
        depthFilterMultiple=3.0,
        snvMaxFilteredBasecallFrac=0.4,
        snvMaxSpanningDeletionFrac=0.75,
        indelMaxWindowFilteredBasecallFrac=0.3,
        ssnvPrior=0.0001,
        sindelPrior=0.000001,
        ssnvNoise=0.0000000005,
        sindelNoiseFactor=2.2,
        ssnvNoiseStrandBiasFrac=0.0,
        minTier1Mapq=20,
        minTier2Mapq=0,
        ssnvQuality_LowerBound=15,
        sindelQuality_LowerBound=40,
        ssnvContamTolerance=0.15,
        indelContamTolerance=0.15,
        isWriteRealignedBam=0,
        docker_image=None
):
    def load_config(file_name):
        config_parser = ConfigParser.ConfigParser()

        config_parser.optionxform = str

        config_parser.read(file_name)

        return dict(config_parser.items('StrelkaSomatic'))

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    region = '{}:{}-{}'.format(*region.split('-'))

    genome_size = sum(known_sizes.values())

    os.makedirs(tmp_dir)

    tmp_indel_file = os.path.join(tmp_dir, 'indels.vcf')

    tmp_snv_file = os.path.join(tmp_dir, 'snvs.vcf')

    stats_file = os.path.join(tmp_dir, 'stats.txt')

    cmd = [
        'run_strelka',
        normal_bam_file,
        tumour_bam_file,
        tmp_indel_file,
        tmp_snv_file,
        stats_file,
        # strelkaSharedWorkflow.py
        region,
        ref_genome_fasta_file,
        genome_size,
        '-max-indel-size', 50,
        # strelkaSomaticWorkflow.py
        '-min-mapping-quality', minTier1Mapq,
        '-min-qscore', 0,
        '-max-window-mismatch', 3, 20,
        '-indel-nonsite-match-prob', 0.5,
        '--somatic-snv-rate', ssnvPrior,
        '--shared-site-error-rate', ssnvNoise,
        '--shared-site-error-strand-bias-fraction', ssnvNoiseStrandBiasFrac,
        '--somatic-indel-rate', sindelPrior,
        '--shared-indel-error-factor', sindelNoiseFactor,
        '--tier2-min-mapping-quality', minTier2Mapq,
        '--tier2-mismatch-density-filter-count', 10,
        '--tier2-indel-nonsite-match-prob', 0.25,
        '--tier2-include-singleton',
        '--tier2-include-anomalous',
        '--strelka-snv-max-filtered-basecall-frac', snvMaxFilteredBasecallFrac,
        '--strelka-snv-max-spanning-deletion-frac', snvMaxSpanningDeletionFrac,
        '--strelka-snv-min-qss-ref', ssnvQuality_LowerBound,
        '--strelka-indel-max-window-filtered-basecall-frac', indelMaxWindowFilteredBasecallFrac,
        '--strelka-indel-min-qsi-ref', sindelQuality_LowerBound,
        '--ssnv-contam-tolerance', ssnvContamTolerance,
        '--indel-contam-tolerance', indelContamTolerance,
    ]

    if not is_exome:
        cmd.extend([
            '--strelka-chrom-depth-file', chrom_depth_file,
            '--strelka-max-depth-factor', depthFilterMultiple,
        ])

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    shutil.move(tmp_indel_file, indel_file)

    shutil.move(tmp_snv_file, snv_file)


def get_known_chromosome_sizes(size_file, chromosomes):
    sizes = {}

    with open(size_file, 'r') as fh:
        reader = csv.DictReader(fh, ['path', 'chrom', 'known_size', 'size'], delimiter='\t')

        for row in reader:
            if row['chrom'] not in chromosomes:
                continue

            sizes[row['chrom']] = int(row['known_size'])

    return sizes


def count_fasta_bases(ref_genome_fasta_file, out_file, docker_image=None):
    cmd = [
        'countFastaBases',
        ref_genome_fasta_file,
        '>',
        out_file
    ]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)
