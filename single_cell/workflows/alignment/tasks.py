import os
import shutil
import warnings
from scripts import CollectMetrics
from scripts import GenerateCNMatrix
from scripts import RunTrimGalore
from scripts import SummaryMetrics

from single_cell.utils import picardutils
from single_cell.utils import bamutils
from single_cell.utils import helpers
from single_cell.utils import csvutils
from single_cell.utils import gatkutils
from single_cell.utils import hdfutils

from single_cell.utils.singlecell_copynumber_plot_utils import PlotMetrics


def plot_metrics(metrics, output, plot_title,gc_matrix, gc_content):
    plot = PlotMetrics(metrics, output, plot_title, gcbias_matrix=gc_matrix, gc_content=gc_content)
    plot.plot_alignment_metrics()


def get_summary_metrics(infile, output):
    summ = SummaryMetrics(infile, output)
    summ.main()


def annotate_metrics(infile, sample_info, outfile):
    
    hdfutils.annotate_store_with_dict(infile, sample_info, outfile)

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
                      config, input_cell_id):
    merge_filenames = []
    for (_, cell_id), filename in input_filenames.iteritems():
        if input_cell_id != cell_id:
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


def get_readgroup(run_id, cell_id, library_id, seqinfo, sample_info):
    platform = 'illumina'
    centre = 'UBCBRC' if seqinfo.lower() == 'nextseq' else 'BCCAGSC'

    barcode = sample_info["i7_index"] + "-" + sample_info["i5_index"]

    read_group_template = (
        '@RG\\tID:' + str(library_id) + '_' + cell_id + '_' + str(run_id) +
        '\\tPL:' + platform +
        '\\tPU:' + str(run_id) +
        '\\tLB:' + str(library_id) + '_' + cell_id +
        '\\tSM:' + cell_id +
        '\\tCN:' + centre +
        '\\tKS:' + barcode)

    return read_group_template



def align_pe(fastq1, fastq2, output, reports, metrics, tempdir,
             reference, source, sample_info, cell_id, lane_id, library_id, config):

    readgroup = get_readgroup(lane_id, cell_id, library_id, source, sample_info)

    run_fastqc(fastq1, fastq2, reports, tempdir)

    aln_temp = os.path.join(tempdir, "temp_alignments.bam")
    if config["aligner"] == "bwa-mem":
        bamutils.bwa_mem_paired_end(fastq1, fastq2, aln_temp, reference, readgroup)
    elif config["aligner"] == "bwa-aln":
        if not source.lower() == "nextseq":
            fastq1, fastq2 = trim_fastqs(fastq1, fastq2, cell_id, tempdir, config)
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
    for cell_id,infile in infiles.iteritems():
        tempout = os.path.join(tempdir, os.path.basename(infile)+".parsed.csv") 
        tempouts.append(tempout)
        gen_gc = GenerateCNMatrix(infile, tempout, ',',
                                  'NORMALIZED_COVERAGE', cell_id,
                                  'gcbias')
        gen_gc.main()

    merged_csv = os.path.join(tempdir, "merged_gc_metrics.csv")
    csvutils.concatenate_csv(tempouts, merged_csv)

    hdfutils.convert_csv_to_hdf(merged_csv, outfile, '/alignment/gc_metrics')

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


    merged_metrics_csv = os.path.join(tempdir, 'merged_alignment_metrics.csv')
    csvutils.concatenate_csv(sample_outputs, merged_metrics_csv)

    hdfutils.convert_csv_to_hdf(merged_metrics_csv, merged_metrics, '/alignment/metrics/')

def run_trimgalore(seq1, seq2, fq_r1, fq_r2, trimgalore, cutadapt, tempdir,
                   adapter, adapter2, report_r1, report_r2, qc_report_r1,
                   qc_report_r2, qc_zip_r1, qc_zip_r2):

    run_tg = RunTrimGalore(seq1, seq2, fq_r1, fq_r2, trimgalore, cutadapt,
                           tempdir, adapter, adapter2, report_r1, report_r2,
                           qc_report_r1, qc_report_r2, qc_zip_r1, qc_zip_r2)
    run_tg.run_trimgalore()
    run_tg.gather_outputs()


def trim_fastqs(fastq1, fastq2, cell_id, tempdir, config):
    """
    run fastqc on both fastq files
    run trimgalore if needed, copy if not.
    """

    trim1 = os.path.join(tempdir, "fastq_R1_trimmed.fastq.gz")
    trim2 = os.path.join(tempdir, "fastq_R2_trimmed.fastq.gz")

    reports_dir = os.path.join(tempdir, 'fastqc_reports')
    if not os.path.exists(reports_dir):
        helpers.makedirs(reports_dir)

    rep1 = os.path.join(reports_dir, '{}_trimgalore_R1.html'.format(cell_id))
    rep2 = os.path.join(reports_dir, '{}_trimgalore_R2.html'.format(cell_id))
    qcrep1 = os.path.join(reports_dir, '{}_trimgalore_qc_R1.html'.format(cell_id))
    qcrep2 = os.path.join(reports_dir, '{}_trimgalore_qc_R2.html'.format(cell_id))
    qczip1 = os.path.join(reports_dir, '{}_trimgalore_qc_R1.zip'.format(cell_id))
    qczip2 = os.path.join(reports_dir, '{}_trimgalore_qc_R2.zip'.format(cell_id))

    run_trimgalore(fastq1, fastq2, trim1, trim2, 'trim_galore', 'cutadapt',
                   tempdir, config['adapter'], config['adapter2'],
                   rep1, rep2, qcrep1, qcrep2, qczip1, qczip2)

    return trim1, trim2

