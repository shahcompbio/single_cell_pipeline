import logging
import os

import pypeliner
import single_cell.workflows.align.fastqscreen as fastqscreen
from single_cell.utils import bamutils
from single_cell.utils import helpers
from single_cell.utils import picardutils

from .scripts import RunTrimGalore


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


def run_fastqc(fastq1, fastq2, reports_dir, tempdir, containers):
    """
    run fastqc on both fastq files
    run trimgalore if needed, copy if not.
    """
    # empty fastq files
    if os.stat(fastq1).st_size < 100 and os.stat(fastq2).st_size < 100:
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

    run_tg = RunTrimGalore(
        fastq1, fastq2, trim1, trim2, 'trim_galore', 'cutadapt', tempdir,
        adapter, adapter2, rep1, rep2, qcrep1, qcrep2, qczip1, qczip2,
        trimgalore_docker
    )
    run_tg.run_trimgalore()
    run_tg.gather_outputs()

    return trim1, trim2


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
        fastq1, fastq2, output, reports_dir, tempdir, reference,
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

    run_fastqc(filtered_fastq_r1, filtered_fastq_r2, reports_dir, tempdir, containers)

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

    metrics = os.path.join(reports_dir, 'flagstat_metrics.txt')
    bamutils.bam_flagstat(output, metrics, docker_image=containers['samtools'])


def extract_mt_chromosome(input_bam, mt_bam, docker_image=None, mt_chrom_name='MT'):
    cmd = ['samtools', 'view', '-bh', input_bam, mt_chrom_name, '>', mt_bam]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    cmd = ['samtools', 'index', mt_bam, mt_bam + '.bai']
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def align_lanes(
        fastq1, fastq2, output, output_mt, reports, tempdir, reference, laneinfo,
        sample_info, cell_id, library_id, aligner, containers, adapter,
        adapter2, fastqscreen_detailed_metrics,
        fastqscreen_summary_metrics, fastqscreen_params, mt_chrom_name='MT'
):
    lane_bams = []
    detailed_counts = []
    summary_counts = []

    for lane_id in fastq1:
        reports_dir = os.path.join(tempdir, 'reports_per_lane', lane_id)

        helpers.makedirs(reports_dir)

        lane_tempdir = os.path.join(tempdir, lane_id, 'lane_temp')
        lane_bam = os.path.join(tempdir, lane_id, 'aligned.bam')

        lane_bams.append(lane_bam)

        screen_detailed = os.path.join(reports_dir, 'detailed.txt')
        screen_summary = os.path.join(reports_dir, 'summary.txt')

        detailed_counts.append(screen_detailed)
        summary_counts.append(screen_summary)

        align_pe(
            fastq1[lane_id], fastq2[lane_id], lane_bam, reports_dir,
            lane_tempdir, reference, laneinfo[lane_id]['trim'],
            laneinfo[lane_id]['center'], sample_info, cell_id, lane_id,
            library_id, aligner, containers, adapter, adapter2,
            screen_detailed, screen_summary, fastqscreen_params
        )

    helpers.make_tarfile(reports, os.path.join(tempdir, 'reports_per_lane'))

    merge_postprocess_bams(lane_bams, output, os.path.join(tempdir, 'merge_bams'), containers)

    fastqscreen.merge_fastq_screen_counts(
        detailed_counts, summary_counts,
        fastqscreen_detailed_metrics, fastqscreen_summary_metrics
    )

    extract_mt_chromosome(output, output_mt, docker_image=containers['samtools'], mt_chrom_name=mt_chrom_name)
