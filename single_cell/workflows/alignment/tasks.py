import pypeliner.commandline
import os
import errno
import shutil
import warnings
import tarfile


def makedirs(directory):
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def produce_fastqc_report(fastq_filename, output_html, output_plots, temp_dir):
    makedirs(temp_dir)

    pypeliner.commandline.execute(
        'fastqc',
        '--outdir=' + temp_dir,
        fastq_filename)

    if fastq_filename.endswith(".fastq.gz"):
        fastq_basename = fastq_filename.replace(".fastq.gz", "")
    if fastq_filename.endswith(".fq.gz"):
        fastq_basename = fastq_filename.replace(".fq.gz", "")
    else:
        raise Exception("Unknown file type")

    output_basename = os.path.join(temp_dir, fastq_basename)

    shutil.move(output_basename + '_fastqc.zip', output_plots)
    shutil.move(output_basename + '_fastqc.html', output_html)


def run_fastqc(fastq1, fastq2, reports, tempdir):
    """
    run fastqc on both fastq files
    run trimgalore if needed, copy if not.
    """
    reports_dir = os.path.join(tempdir, 'fastqc_reports')
    if not os.path.exists(reports_dir):
        makedirs(reports_dir)

    out_html = os.path.join(reports_dir, 'fastqc_R1.html')
    out_plot = os.path.join(reports_dir, 'fastqc_R1.zip')
    if not os.path.getsize(fastq1) == 0:
        produce_fastqc_report(fastq1, out_html, out_plot, tempdir)
    else:
        warnings.warn("fastq file %s is empty, skipping fastqc" % fastq1)

    out_html = os.path.join(reports_dir, 'fastqc_R2.html')
    out_plot = os.path.join(reports_dir, 'fastqc_R2.zip')
    if not os.path.getsize(fastq2) == 0:
        produce_fastqc_report(fastq2, out_html, out_plot, tempdir)
    else:
        warnings.warn("fastq file %s is empty, skipping fastqc" % fastq1)

    make_tarfile(reports, reports_dir)


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


def bam_sort(bam_filename, sorted_bam_filename):
    pypeliner.commandline.execute(
        'picard', '-Xmx2G', '-Xms2G',
        '-XX:ParallelGCThreads=1',
        'SortSam',
        'INPUT=' + bam_filename,
        'OUTPUT=' + sorted_bam_filename,
        'SORT_ORDER=coordinate',
        'VALIDATION_STRINGENCY=LENIENT',
        'MAX_RECORDS_IN_RAM=150000')


def bwa_align_paired_end(fastq1, fastq2, output,
                         reference, readgroup):
    """
    run bwa aln on both fastq files,
    bwa sampe to align, and convert to bam with samtools view
    """
    pypeliner.commandline.execute(
        'bwa', 'mem', '-M', '-R', readgroup,
        reference, fastq1, fastq2, '|',
        'samtools', 'view', '-bSh', '-',
        '>', output,
    )


def run_flagstat(bam, metrics):

    pypeliner.commandline.execute(
        'samtools', 'flagstat',
        bam,
        '>',
        metrics
    )


def align_pe(fastq1, fastq2, output, reports, metrics, tempdir,
             reference, source, sample_id, lane_id, library_id):

    readgroup = get_readgroup(lane_id, sample_id, library_id, source)

    run_fastqc(fastq1, fastq2, reports, tempdir)

    aln_temp = os.path.join(tempdir, "temp_alignments.bam")
    bwa_align_paired_end(fastq1, fastq2, aln_temp, reference, readgroup)

    bam_sort(aln_temp, output)

    run_flagstat(output, metrics)
