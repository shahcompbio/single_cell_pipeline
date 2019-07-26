import logging
import os
import shutil

from single_cell.utils import bamutils
from single_cell.utils import csvutils
from single_cell.utils import gatkutils
from single_cell.utils import helpers
from single_cell.utils import picardutils
from single_cell.utils.singlecell_copynumber_plot_utils import PlotMetrics

from scripts import CollectMetrics
from scripts import GenerateCNMatrix
from scripts import RunTrimGalore
from scripts import SummaryMetrics


def merge_bams(inputs, output, output_index, containers):
    picardutils.merge_bams(inputs, output, docker_image=containers['picard'])
    bamutils.bam_index(output, output_index, docker_image=containers['samtools'])


def merge_realignment(input_filenames, output_filename,
                      config, input_cell_id):
    merge_filenames = []
    for (_, cell_id), filename in input_filenames.iteritems():
        if input_cell_id != cell_id:
            continue
        merge_filenames.append(filename)

    merge_bams(merge_filenames, output_filename, output_filename + ".bai")


def realign(input_bams, input_bais, output_bams, tempdir, config, interval):
    container_ctx = helpers.get_container_ctx(config['containers'], 'samtools', docker_only=True)

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
    gatkutils.generate_targets(input_bams, config, targets, interval, **container_ctx)

    # run gatk realigner
    gatkutils.gatk_realigner(new_inputs, config, targets, interval, tempdir, **container_ctx)

    # copy generated files in temp dir to the specified output paths
    for key in input_bams.keys():
        realigned_bam = os.path.join(tempdir, key + '_indel_realigned.bam')
        realigned_bai = os.path.join(tempdir, key + '_indel_realigned.bai')
        output_bam_filename = output_bams[key]
        output_bai_filename = output_bam_filename + '.bai'

        shutil.move(realigned_bam, output_bam_filename)
        shutil.move(realigned_bai, output_bai_filename)


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
        containers, adapter, adapter2, fastqscreen_params
):
    readgroup = get_readgroup(
        lane_id, cell_id, library_id, centre, sample_info
    )

    run_fastqc(fastq1, fastq2, reports, tempdir, containers)

    aln_temp = os.path.join(tempdir, "temp_alignments.bam")

    if aligner == "bwa-aln" and trim:
        fastq1, fastq2 = trim_fastqs(
            fastq1, fastq2, cell_id, tempdir,
            adapter, adapter2, containers['trimgalore']
        )

    align_pe_with_bwa(
        fastq1, fastq2, aln_temp, reference, readgroup,
        tempdir, containers, aligner=aligner
    )

    picardutils.bam_sort(aln_temp, output, tempdir, docker_image=containers['picard'])

    bamutils.bam_flagstat(output, metrics, docker_image=containers['samtools'])


def postprocess_bam(infile, outfile, tempdir, containers):
    outfile_index = outfile + '.bai'

    if not os.path.exists(tempdir):
        helpers.makedirs(tempdir)

    sorted_bam = os.path.join(tempdir, 'sorted.bam')
    picardutils.bam_sort(infile, sorted_bam, tempdir,
                         docker_image=containers['picard'])

    markdups_metrics = os.path.join(tempdir, 'markdups_metrics.txt')
    picardutils.bam_markdups(sorted_bam, outfile, markdups_metrics, tempdir,
                             docker_image=containers['picard'])

    bamutils.bam_index(outfile, outfile_index, docker_image=containers['samtools'])


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
    for cell_id, infile in infiles.iteritems():
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
