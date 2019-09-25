from __future__ import division

import logging
import os

import single_cell.workflows.align.fastqscreen as fastqscreen
from single_cell.utils import bamutils
from single_cell.utils import csvutils
from single_cell.utils import helpers
from single_cell.utils import picardutils
from single_cell.utils.singlecell_copynumber_plot_utils import PlotMetrics

from .scripts import CollectMetrics
from .scripts import GenerateCNMatrix
from .scripts import RunTrimGalore
from .scripts import SummaryMetrics


class LibraryContaminationError(Exception):
    pass


def add_contamination_status(
        infile, outfile,
        reference='grch37', ref_threshold=0.6, alt_threshold=0.2,
        strict_validation=True
):
    data = csvutils.read_csv_and_yaml(infile)

    data = data.set_index('cell_id', drop=False)

    fastqscreen_cols = [col for col in data.columns.values if col.startswith('fastqscreen_')]

    reference = "fastqscreen_{}".format(reference)
    if reference not in fastqscreen_cols:
        raise Exception("Could not find the fastq screen counts")

    alts = [col for col in fastqscreen_cols if not col == reference]

    data['is_contaminated'] = False

    perc_ref = data[reference] / data['total_reads']
    data.loc[perc_ref <= ref_threshold, 'is_contaminated'] = True

    for altcol in alts:
        perc_alt = data[altcol] / data['total_reads']
        data.loc[perc_alt > alt_threshold, 'is_contaminated'] = True

    csvutils.write_dataframe_to_csv_and_yaml(
        data, outfile, write_header=True
    )

    # get cells that are contaminated and have enopugh human reads
    check_df = data.loc[data['is_contaminated'] == True]
    check_df['perc_ref'] = data[reference] / data['total_reads']
    check_df = check_df[check_df['perc_ref'] > ref_threshold]
    if strict_validation and (len(check_df) / len(data) > 0.2):
        raise LibraryContaminationError("over 20% of cells are contaminated")


def merge_postprocess_bams(inputs, output, tempdir, containers):
    helpers.makedirs(tempdir)
    merged_out = os.path.join(tempdir, 'merged_lanes.bam')

    picardutils.merge_bams(inputs, merged_out, docker_image=containers['picard'])
    bamutils.bam_index(merged_out, merged_out + '.bai', docker_image=containers['samtools'])

    sorted_bam = os.path.join(tempdir, 'sorted.bam')
    picardutils.bam_sort(merged_out, sorted_bam, tempdir,
                         docker_image=containers['picard'])

    markdups_metrics = os.path.join(tempdir, 'markdups_metrics.txt')
    picardutils.bam_markdups(sorted_bam, output, markdups_metrics, tempdir,
                             docker_image=containers['picard'])

    bamutils.bam_index(output, output + '.bai', docker_image=containers['samtools'])


def run_fastqc(fastq1, fastq2, reports, tempdir, containers):
    """
    run fastqc on both fastq files
    run trimgalore if needed, copy if not.
    """
    reports_dir = os.path.join(tempdir, 'fastqc_reports')
    if not os.path.exists(reports_dir):
        helpers.makedirs(reports_dir)

    # empty fastq files
    if os.stat(fastq1).st_size < 100 and os.stat(fastq2).st_size < 100:
        helpers.make_tarfile(reports, reports_dir)
        return

    out_html = os.path.join(reports_dir, 'fastqc_R1.html')
    out_plot = os.path.join(reports_dir, 'fastqc_R1.zip')
    if not os.path.getsize(fastq1) == 0:
        bamutils.produce_fastqc_report(fastq1, out_html, out_plot, tempdir,
                                       docker_image=containers['fastqc'])
    else:
        logging.getLogger("single_cell.align.tasks").warn(
            "fastq file %s is empty, skipping fastqc" % fastq1)

    out_html = os.path.join(reports_dir, 'fastqc_R2.html')
    out_plot = os.path.join(reports_dir, 'fastqc_R2.zip')
    if not os.path.getsize(fastq2) == 0:
        bamutils.produce_fastqc_report(fastq2, out_html, out_plot, tempdir,
                                       docker_image=containers['fastqc'])
    else:
        logging.getLogger("single_cell.align.tasks").warn(
            "fastq file %s is empty, skipping fastqc" % fastq1)

    helpers.make_tarfile(reports, reports_dir)


def get_readgroup(run_id, cell_id, library_id, centre, sample_info):
    platform = 'illumina'

    if not centre:
        logging.getLogger("single_cell.align.tasks").warn("no sequencing centre specified")
        centre = "NA"

    barcode = sample_info["primer_i7"] + "-" + sample_info["primer_i5"]

    read_group_template = (
            r'@RG\tID:' + str(library_id) + '_' + cell_id + '_' + str(run_id) +
            r'\tPL:' + platform +
            r'\tPU:' + str(run_id) +
            r'\tLB:' + str(library_id) + '_' + cell_id +
            r'\tSM:' + cell_id +
            r'\tCN:' + centre +
            r'\tKS:' + barcode)

    return read_group_template


def align_pe_with_bwa(
        fastq1, fastq2, output, reference, readgroup, tempdir,
        containers, aligner='bwa-aln'
):
    samfile = os.path.join(tempdir, "bwamem.sam")

    if aligner == 'bwa-aln':
        bamutils.bwa_aln_paired_end(fastq1, fastq2, samfile, tempdir, reference, readgroup,
                                    docker_image=containers['bwa'])
    elif aligner == 'bwa-mem':
        bamutils.bwa_mem_paired_end(fastq1, fastq2, samfile, reference, readgroup,
                                    docker_image=containers['bwa'])
    else:
        raise Exception(
            "Aligner %s not supported, pipeline supports bwa-aln and bwa-mem" %
            aligner)

    bamutils.samtools_sam_to_bam(samfile, output,
                                 docker_image=containers['samtools'])


def align_pe(
        fastq1, fastq2, output, reports, metrics, tempdir, reference,
        trim, centre, sample_info, cell_id, lane_id, library_id, aligner,
        containers, adapter, adapter2,
        fastqscreen_detailed_metrics, fastqscreen_summary_metrics,
        fastqscreen_params
):
    fastqscreen_tempdir = os.path.join(tempdir, 'fastq_screen')
    helpers.makedirs(fastqscreen_tempdir)

    filtered_fastq_r1 = os.path.join(fastqscreen_tempdir, "fastq_r1.fastq.gz")
    filtered_fastq_r2 = os.path.join(fastqscreen_tempdir, "fastq_r2.fastq.gz")

    fastqscreen.organism_filter(
        fastq1, fastq2, filtered_fastq_r1, filtered_fastq_r2,
        fastqscreen_detailed_metrics, fastqscreen_summary_metrics,
        fastqscreen_tempdir, cell_id, fastqscreen_params,
        reference, docker_image=containers['fastq_screen'],
        filter_contaminated_reads=fastqscreen_params['filter_contaminated_reads'],
    )

    readgroup = get_readgroup(
        lane_id, cell_id, library_id, centre, sample_info
    )

    run_fastqc(filtered_fastq_r1, filtered_fastq_r2, reports, tempdir, containers)

    aln_temp = os.path.join(tempdir, "temp_alignments.bam")

    if aligner == "bwa-aln" and trim:
        filtered_fastq_r1, filtered_fastq_r2 = trim_fastqs(
            filtered_fastq_r1, filtered_fastq_r2, cell_id, tempdir,
            adapter, adapter2, containers['trimgalore']
        )

    align_pe_with_bwa(
        filtered_fastq_r1, filtered_fastq_r2, aln_temp, reference, readgroup,
        tempdir, containers, aligner=aligner
    )

    picardutils.bam_sort(aln_temp, output, tempdir, docker_image=containers['picard'])

    bamutils.bam_flagstat(output, metrics, docker_image=containers['samtools'])


def run_trimgalore(seq1, seq2, fq_r1, fq_r2, trimgalore, cutadapt, tempdir,
                   adapter, adapter2, report_r1, report_r2, qc_report_r1,
                   qc_report_r2, qc_zip_r1, qc_zip_r2, docker_image):
    run_tg = RunTrimGalore(seq1, seq2, fq_r1, fq_r2, trimgalore, cutadapt,
                           tempdir, adapter, adapter2, report_r1, report_r2,
                           qc_report_r1, qc_report_r2, qc_zip_r1, qc_zip_r2,
                           docker_image)
    run_tg.run_trimgalore()
    run_tg.gather_outputs()


def trim_fastqs(fastq1, fastq2, cell_id, tempdir, adapter, adapter2, trimgalore_docker):
    """
    run fastqc on both fastq files
    run trimgalore if needed, copy if not.
    """
    with helpers.getFileHandle(fastq1) as reader:
        if not reader.readline():
            return fastq1, fastq2

    trim1 = os.path.join(tempdir, "fastq_R1_trimmed.fastq.gz")
    trim2 = os.path.join(tempdir, "fastq_R2_trimmed.fastq.gz")

    reports_dir = os.path.join(tempdir, 'fastqc_reports')
    if not os.path.exists(reports_dir):
        helpers.makedirs(reports_dir)

    rep1 = os.path.join(reports_dir, '{}_trimgalore_R1.html'.format(cell_id))
    rep2 = os.path.join(reports_dir, '{}_trimgalore_R2.html'.format(cell_id))
    qcrep1 = os.path.join(
        reports_dir,
        '{}_trimgalore_qc_R1.html'.format(cell_id))
    qcrep2 = os.path.join(
        reports_dir,
        '{}_trimgalore_qc_R2.html'.format(cell_id))
    qczip1 = os.path.join(
        reports_dir,
        '{}_trimgalore_qc_R1.zip'.format(cell_id))
    qczip2 = os.path.join(
        reports_dir,
        '{}_trimgalore_qc_R2.zip'.format(cell_id))

    run_trimgalore(fastq1, fastq2, trim1, trim2, 'trim_galore', 'cutadapt',
                   tempdir, adapter, adapter2,
                   rep1, rep2, qcrep1, qcrep2, qczip1, qczip2, trimgalore_docker)

    return trim1, trim2


def plot_metrics(metrics, output, plot_title, gc_matrix, gc_content):
    plot = PlotMetrics(
        metrics, output, plot_title, gcbias_matrix=gc_matrix,
        gc_content=gc_content, tablename='/alignment/metrics',
        gc_tablename='/alignment/gc_metrics')
    plot.plot_alignment_metrics()


def get_summary_metrics(infile, output):
    summ = SummaryMetrics(infile, output)
    summ.main()


def collect_gc(infiles, outfile, tempdir):
    helpers.makedirs(tempdir)

    tempouts = []
    for cell_id, infile in infiles.items():
        tempout = os.path.join(
            tempdir,
            "{}.parsed.csv".format(cell_id))
        tempouts.append(tempout)
        gen_gc = GenerateCNMatrix(infile, tempout, ',',
                                  'NORMALIZED_COVERAGE', cell_id,
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
        outfile = os.path.join(tempdir, sample + "_metrics.csv")
        sample_outputs.append(outfile)

        collmet = CollectMetrics(wgs, insrt, flgstat,
                                 mkdup, outfile, sample)
        collmet.main()

    csvutils.concatenate_csv(sample_outputs, merged_metrics)


def picard_wgs_dup(
        input_bam, markdups_bam, markdups_metrics, tempdir,
        ref_genome, wgs_metrics, picard_wgs_params,
        picard_docker=None

):
    tempdir_markdups = os.path.join(tempdir, 'markdups')
    helpers.makedirs(tempdir_markdups)

    picardutils.bam_markdups(
        input_bam,
        markdups_bam,
        markdups_metrics,
        tempdir_markdups,
        docker_image=picard_docker
    )

    tempdir_wgs = os.path.join(tempdir, 'wgs')
    helpers.makedirs(tempdir_wgs)

    picardutils.bam_collect_wgs_metrics(
        input_bam,
        ref_genome,
        wgs_metrics,
        picard_wgs_params,
        tempdir_wgs,
        docker_image=picard_docker
    )


def picard_insert_gc_flagstat(
        input_bam, ref_genome, gc_metrics, gc_metrics_summary, gc_metrics_pdf,
        tempdir, flagstat_metrics, insert_metrics, insert_pdf, picard_docker=None,
        samtools_docker=None
):
    bamutils.bam_flagstat(
        input_bam,
        flagstat_metrics,
        docker_image=samtools_docker
    )

    gc_tempdir = os.path.join(tempdir, 'gc')
    helpers.makedirs(gc_tempdir)

    picardutils.bam_collect_gc_metrics(
        input_bam,
        ref_genome,
        gc_metrics,
        gc_metrics_summary,
        gc_metrics_pdf,
        gc_tempdir,
        docker_image=picard_docker
    )

    insert_tempdir = os.path.join(tempdir, 'insert')
    helpers.makedirs(insert_tempdir)
    picardutils.bam_collect_insert_metrics(
        input_bam,
        flagstat_metrics,
        insert_metrics,
        insert_pdf,
        insert_tempdir,
        docker_image=picard_docker
    )
