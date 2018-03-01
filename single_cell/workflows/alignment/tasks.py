import pypeliner.commandline
import os
import errno
import shutil
import warnings
import tarfile
import pandas as pd
from scripts import CollectMetrics
from scripts import GenerateCNMatrix

from single_cell.utils import picardutils
from single_cell.utils import bamutils
from single_cell.utils import helpers
from single_cell.utils import csvutils
from single_cell.utils import gatkutils


def merge_all_metrics(infiles, outfile):
    csvutils.merge_csv(infiles, outfile, "outer", "cell_id")

def bam_collect_wgs_metrics(bam_filename, ref_genome, metrics_filename, config, tempdir):
    picardutils.bam_collect_wgs_metrics(bam_filename, ref_genome, metrics_filename, config, tempdir)

def bam_collect_gc_metrics(bam_filename, ref_genome, metrics_filename, summary_filename, chart_filename, tempdir):
    picardutils.bam_collect_gc_metrics(bam_filename, ref_genome, metrics_filename, summary_filename, chart_filename, tempdir)

def bam_collect_insert_metrics(bam_filename, flagstat_metrics_filename, metrics_filename, histogram_filename, tempdir):
    picardutils.bam_collect_insert_metrics(bam_filename, flagstat_metrics_filename, metrics_filename, histogram_filename, tempdir)

def merge_bams(inputs, output, output_index):

    picardutils.merge_bams(inputs, output)
    bamutils.bam_index(output, output_index)

def merge_realignment(input_filenames, output_filename,
                      config, input_sample_id):
    merge_filenames = []
    for (_, sample_id), filename in input_filenames.iteritems():
        if input_sample_id != sample_id:
            continue
        merge_filenames.append(filename)

    merge_bams(merge_filenames, output_filename, output_filename+".bai")

def realign(input_bams, input_bais, output_bams, tempdir, config, interval):

    # make the dir
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    # symlink inputs to tempdir, inputs have same filename but they should be
    # different for mapping file nwayout to work
    # realign
    new_inputs = {}
    for key, bamfile in input_bams.iteritems():
        new_bam = os.path.join(tempdir, key + '.bam')
        new_bai = os.path.join(tempdir, key + '.bam.bai')

        shutil.copy(bamfile, new_bam)
        shutil.copy(bamfile + '.bai', new_bai)
        new_inputs[key] = new_bam

    # save intervals file in tempdir
    targets = os.path.join(tempdir, 'realn_positions.intervals')
    gatkutils.generate_targets(input_bams, config, targets, interval)

    # run gatk realigner
    gatkutils.gatk_realigner(new_inputs, config, targets, interval, tempdir)

    # copy generated files in temp dir to the specified output paths
    for key in input_bams.keys():
        realigned_bam = os.path.join(tempdir, key + '_indel_realigned.bam')
        realigned_bai = os.path.join(tempdir, key + '_indel_realigned.bai')
        output_bam_filename = output_bams[key]
        output_bai_filename = output_bam_filename + '.bai'

        shutil.move(realigned_bam, output_bam_filename)
        shutil.move(realigned_bai, output_bai_filename)


def run_fastqc(fastq1, fastq2, reports, tempdir):
    """
    run fastqc on both fastq files
    run trimgalore if needed, copy if not.
    """
    reports_dir = os.path.join(tempdir, 'fastqc_reports')
    if not os.path.exists(reports_dir):
        helpers.makedirs(reports_dir)

    out_html = os.path.join(reports_dir, 'fastqc_R1.html')
    out_plot = os.path.join(reports_dir, 'fastqc_R1.zip')
    if not os.path.getsize(fastq1) == 0:
        bamutils.produce_fastqc_report(fastq1, out_html, out_plot, tempdir)
    else:
        warnings.warn("fastq file %s is empty, skipping fastqc" % fastq1)

    out_html = os.path.join(reports_dir, 'fastqc_R2.html')
    out_plot = os.path.join(reports_dir, 'fastqc_R2.zip')
    if not os.path.getsize(fastq2) == 0:
        bamutils.produce_fastqc_report(fastq2, out_html, out_plot, tempdir)
    else:
        warnings.warn("fastq file %s is empty, skipping fastqc" % fastq1)

    helpers.make_tarfile(reports, reports_dir)


def get_readgroup(run_id, sample_id, library_id, seqinfo):
    platform = 'illumina'
    centre = 'UBCBRC' if seqinfo == 'nextseq' else 'BCCAGSC'

    read_group_template = (
        '@RG\\tID:' + str(library_id) + '_' + sample_id + '_' + str(run_id) +
        '\\tPL:' + platform +
        '\\tPU:' + str(run_id) +
        '\\tLB:' + str(library_id) + '_' + sample_id +
        '\\tSM:' + sample_id +
        '\\tCN:' + centre)

    return read_group_template



def align_pe(fastq1, fastq2, output, reports, metrics, tempdir,
             reference, source, sample_id, lane_id, library_id, config):

    readgroup = get_readgroup(lane_id, sample_id, library_id, source)

    run_fastqc(fastq1, fastq2, reports, tempdir)

    aln_temp = os.path.join(tempdir, "temp_alignments.bam")
    if config["aligner"] == "bwa-mem":
        bamutils.bwa_mem_paired_end(fastq1, fastq2, aln_temp, reference, readgroup)
    elif config["aligner"] == "bwa-aln":
        bamutils.bwa_aln_paired_end(fastq1, fastq2, aln_temp, tempdir, reference, readgroup)
    else:
        raise Exception("Aligner %s not supported, pipeline supports bwa-aln and bwa-mem" %config["aligner"])

    picardutils.bam_sort(aln_temp, output, tempdir)

    bamutils.bam_flagstat(output, metrics)


def postprocess_bam(infile, outfile, outfile_index, tempdir,
                    config, markdups_metrics, flagstat_metrics):
    
    if not os.path.exists(tempdir):
        helpers.makedirs(tempdir)
    
    sorted_bam = os.path.join(tempdir, 'sorted.bam')
    
    picardutils.bam_sort(infile, sorted_bam, tempdir)

    picardutils.bam_markdups(sorted_bam, outfile, markdups_metrics, tempdir)

    bamutils.bam_index(outfile, outfile_index)
    
    bamutils.bam_flagstat(outfile, flagstat_metrics)


def collect_gc(infiles, outfile, tempdir):

    helpers.makedirs(tempdir)

    tempouts = []
    for sample_id,infile in infiles.iteritems():
        tempout = os.path.join(tempdir, os.path.basename(infile)+".parsed.csv") 
        tempouts.append(tempout)
        gen_gc = GenerateCNMatrix(infile, tempout, ',',
                                  'NORMALIZED_COVERAGE', sample_id,
                                  'gcbias')
        gen_gc.main()

    csvutils.concatenate_csv(tempouts, outfile)

def collect_metrics(flagstat_metrics, markdups_metrics, insert_metrics,
                    wgs_metrics, tempdir, merged_metrics):

    helpers.makedirs(tempdir)
    sample_outputs = []
    for sample in flagstat_metrics.keys():
        flgstat = flagstat_metrics[sample]
        mkdup = markdups_metrics[sample]
        insrt = insert_metrics[sample]
        wgs = wgs_metrics[sample]
        outfile = os.path.join(tempdir, sample+"_metrics.csv")
        sample_outputs.append(outfile)

        collmet = CollectMetrics(wgs, insrt, flgstat,
                                 mkdup, outfile, sample)
        collmet.main()

    csvutils.concatenate_csv(sample_outputs, merged_metrics)

