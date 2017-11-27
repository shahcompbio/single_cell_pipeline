import pypeliner.commandline
import os
import errno
import shutil
import time
import warnings
from scripts import RunTrimGalore
import tarfile

def copy_files(inputs, outputs):
    for inp, outp in zip(inputs, outputs):
        shutil.copy(inp, outp)


def makedirs(directory):
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def produce_fastqc_report(
        fastq_filename, output_html, output_plots, temp_directory, fastqc):
    makedirs(temp_directory)

    pypeliner.commandline.execute(
        fastqc,
        '--outdir=' + temp_directory,
        fastq_filename)

    fastq_basename = os.path.basename(fastq_filename).split('.')[0]
    output_basename = os.path.join(temp_directory, fastq_basename)

    shutil.move(output_basename + '_fastqc.zip', output_plots)
    shutil.move(output_basename + '_fastqc.html', output_html)


def run_trimgalore(seq1, seq2, fq_r1, fq_r2, trimgalore, cutadapt, tempdir,
                   adapter, adapter2, report_r1, report_r2, qc_report_r1,
                   qc_report_r2, qc_zip_r1, qc_zip_r2):

    run_tg = RunTrimGalore(seq1, seq2, fq_r1, fq_r2, trimgalore, cutadapt,
                           tempdir, adapter, adapter2, report_r1, report_r2,
                           qc_report_r1, qc_report_r2, qc_zip_r1, qc_zip_r2)
    run_tg.run_trimgalore()
    run_tg.gather_outputs()



def trim_fastqs(fastq1, fastq2, trim1, trim2, reports, sample_id, tempdir, source, config):
    """
    run fastqc on both fastq files
    run trimgalore if needed, copy if not.
    """


    reports_dir = os.path.join(tempdir, 'fastqc_reports')
    if not os.path.exists(reports_dir):
        makedirs(reports_dir)

    
    out_html = os.path.join(reports_dir, '{}_fastqc_R1.html'.format(sample_id))
    out_plot = os.path.join(reports_dir, '{}_fastqc_R1.zip'.format(sample_id))

    if not os.path.getsize(fastq1) == 0:
        warnings.warn("fastq file %s is empty, skipping fastqc" %fastq1)
        produce_fastqc_report(fastq1, out_html, out_plot, tempdir, 'fastqc')

    out_html = os.path.join(reports_dir, '{}_fastqc_R2.html'.format(sample_id))
    out_plot = os.path.join(reports_dir, '{}_fastqc_R2.zip'.format(sample_id))

    if not os.path.getsize(fastq2) == 0:
        warnings.warn("fastq file %s is empty, skipping fastqc" %fastq2)
        produce_fastqc_report(fastq2, out_html, out_plot, tempdir, 'fastqc')

    if source == 'hiseq':
        rep1 = os.path.join(reports_dir, '{}_trimgalore_R1.html'.format(sample_id))
        rep2 = os.path.join(reports_dir, '{}_trimgalore_R2.html'.format(sample_id))
        qcrep1 = os.path.join(reports_dir, '{}_trimgalore_qc_R1.html'.format(sample_id))
        qcrep2 = os.path.join(reports_dir, '{}_trimgalore_qc_R2.html'.format(sample_id))
        qczip1 = os.path.join(reports_dir, '{}_trimgalore_qc_R1.zip'.format(sample_id))
        qczip2 = os.path.join(reports_dir, '{}_trimgalore_qc_R2.zip'.format(sample_id))

        
        run_trimgalore(fastq1, fastq2, trim1, trim2, 'trim_galore', 'cutadapt',
                       tempdir, config['adapter'], config['adapter2'],
                       rep1, rep2, qcrep1, qcrep2, qczip1, qczip2)
    else:
        copy_files([fastq1, fastq2], [trim1, trim2])

    make_tarfile(reports, reports_dir)




def get_readgroup(run_id, sample_id, args, config, seqinfo):
    platform = 'illumina'
    centre = 'UBCBRC' if seqinfo[sample_id] == 'nextseq' else 'BCCAGSC'

    if 'read_group' in config:
        if config['read_group']['PL']:
            platform = str(config['read_group']['PL'])
        if config['read_group']['CN']:
            centre = str(config['read_group']['CN'])

    library_id = args['library_id']
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


def bwa_align_paired_end(fastq1, fastq2, output, tempdir,
                     reference, readgroup):
    """
    run bwa aln on both fastq files,
    bwa sampe to align, and convert to bam with samtools view
    """
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    read_1_sai = os.path.join(tempdir, 'read_1.sai')
    read_2_sai = os.path.join(tempdir, 'read_2.sai')

    pypeliner.commandline.execute(
        'bwa',
        'aln',
        reference,
        fastq1,
        '>',
        read_1_sai
    )

    pypeliner.commandline.execute(
        'bwa',
        'aln',
        reference,
        fastq2,
        '>',
        read_2_sai,

    )

    pypeliner.commandline.execute(
        'bwa', 'sampe',
        '-r', readgroup,
        reference,
        read_1_sai,
        read_2_sai,
        fastq1,
        fastq2,
        '|',
        'samtools', 'view',
        '-bSh', '-',
        '>',
        output,
    )

def run_flagstat(bam, metrics):
    
    pypeliner.commandline.execute(
            'samtools', 'flagstat',
            bam,
            '>',
            metrics
        )


def align_pe(fastq1, fastq2, trim1, trim2, output, reports, metrics, tempdir,
             reference, config, readgroup, source, sample_id):


    trim_fastqs(fastq1, fastq2, trim1, trim2, reports, sample_id, tempdir, source, config)

    aln_temp = os.path.join(tempdir, "temp_alignments.bam")
    bwa_align_paired_end(trim1, trim2, aln_temp, tempdir, reference,
                         readgroup)

    time.sleep(30)

    bam_sort(aln_temp, output)

    run_flagstat(output, metrics)

