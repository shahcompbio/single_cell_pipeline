'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import errno
import pypeliner
from scripts import CollectMetrics

from utils import bamutils
from utils import picardutils

def bam_collect_wgs_metrics(bam_filename, ref_genome, metrics_filename, config, tempdir):
    picardutils.bam_collect_wgs_metrics(bam_filename, ref_genome, metrics_filename, config, tempdir, mem="12G")

def bam_collect_gc_metrics(bam_filename, ref_genome, metrics_filename, summary_filename, chart_filename, tempdir):
    picardutils.bam_collect_gc_metrics(bam_filename, ref_genome, metrics_filename, summary_filename, chart_filename, tempdir, mem="12G")

def bam_collect_insert_metrics(bam_filename, flagstat_metrics_filename, metrics_filename, histogram_filename, tempdir):
    picardutils.bam_collect_insert_metrics(bam_filename, flagstat_metrics_filename, metrics_filename, histogram_filename, tempdir, mem="12G")

def merge_bams(inputs, output, output_index):
    picardutils.merge_bams(inputs, output, mem="12G")
    bamutils.bam_index(output, output_index)

def bam_sort(bam_filename, sorted_bam_filename, tmp_dir):
    picardutils.bam_sort(bam_filename, sorted_bam_filename, tmp_dir, mem="12G")

def bam_markdups(bam_filename, markduped_bam_filename, metrics_filename, tmp_dir):
    picardutils.bam_markdups(bam_filename, markduped_bam_filename, metrics_filename, tmp_dir, mem="12G")

def index_bam(infile, outfile):
    bamutils.bam_index(infile, outfile)

def flagstat_bam(infile, outfile):
    bamutils.bam_flagstat(infile, outfile)

def collect_metrics(flagstat_metrics, markdups_metrics, insert_metrics,
                    wgs_metrics, samplesheet, output, sample_id):

    collmet = CollectMetrics(wgs_metrics, insert_metrics, flagstat_metrics,
                             markdups_metrics, output, sample_id)
    collmet.main()
