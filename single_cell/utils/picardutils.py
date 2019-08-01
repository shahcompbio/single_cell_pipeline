'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os

import pypeliner.commandline

from single_cell.utils.helpers import makedirs


def merge_bams(inputs, output, mem="2G", **kwargs):
    if isinstance(inputs, dict):
        inputs = inputs.values()

    cmd = ['picard', '-Xmx' + mem, '-Xms' + mem,
           '-XX:ParallelGCThreads=1',
           'MergeSamFiles',
           'OUTPUT=' + output,
           'SORT_ORDER=coordinate',
           'ASSUME_SORTED=true',
           'VALIDATION_STRINGENCY=LENIENT',
           'MAX_RECORDS_IN_RAM=150000'
           ]

    for bamfile in inputs:
        cmd.append('I=' + os.path.abspath(bamfile))

    pypeliner.commandline.execute(*cmd, **kwargs)


def bam_sort(bam_filename, sorted_bam_filename, tempdir, mem="2G", docker_image=None):
    if not os.path.exists(tempdir):
        makedirs(tempdir)

    pypeliner.commandline.execute(
        'picard', '-Xmx' + mem, '-Xms' + mem,
        '-XX:ParallelGCThreads=1',
        'SortSam',
                  'INPUT=' + bam_filename,
                  'OUTPUT=' + sorted_bam_filename,
        'SORT_ORDER=coordinate',
        'VALIDATION_STRINGENCY=LENIENT',
                  'TMP_DIR=' + tempdir,
        'MAX_RECORDS_IN_RAM=150000',
        docker_image=docker_image)


def bam_markdups(bam_filename, markduped_bam_filename, metrics_filename,
                 tempdir, mem="2G", docker_image=None):
    if not os.path.exists(tempdir):
        makedirs(tempdir)

    pypeliner.commandline.execute(
        'picard', '-Xmx' + mem, '-Xms' + mem,
        '-XX:ParallelGCThreads=1',
        'MarkDuplicates',
                  'INPUT=' + bam_filename,
                  'OUTPUT=' + markduped_bam_filename,
                  'METRICS_FILE=' + metrics_filename,
        'REMOVE_DUPLICATES=False',
        'ASSUME_SORTED=True',
        'VALIDATION_STRINGENCY=LENIENT',
                  'TMP_DIR=' + tempdir,
        'MAX_RECORDS_IN_RAM=150000',
        docker_image=docker_image
    )


def bam_collect_wgs_metrics(bam_filename, ref_genome, metrics_filename,
                            config, tempdir, mem="2G", docker_image=None):
    if not os.path.exists(tempdir):
        makedirs(tempdir)

    pypeliner.commandline.execute(
        'picard', '-Xmx' + mem, '-Xms' + mem,
        '-XX:ParallelGCThreads=1',
        'CollectWgsMetrics',
                  'INPUT=' + bam_filename,
                  'OUTPUT=' + metrics_filename,
                  'REFERENCE_SEQUENCE=' + ref_genome,
                  'MINIMUM_BASE_QUALITY=' +
                  str(config['min_bqual']),
                  'MINIMUM_MAPPING_QUALITY=' +
                  str(config['min_mqual']),
        'COVERAGE_CAP=500',
        'VALIDATION_STRINGENCY=LENIENT',
                  'COUNT_UNPAIRED=' +
                  ('True' if config['count_unpaired'] else 'False'),
                  'TMP_DIR=' + tempdir,
        'MAX_RECORDS_IN_RAM=150000',
        docker_image=docker_image
    )


def bam_collect_gc_metrics(bam_filename, ref_genome, metrics_filename,
                           summary_filename, chart_filename, tempdir,
                           mem="2G", docker_image=None):
    if not os.path.exists(tempdir):
        makedirs(tempdir)

    pypeliner.commandline.execute(
        'picard', '-Xmx' + mem, '-Xms' + mem,
        '-XX:ParallelGCThreads=1',
        'CollectGcBiasMetrics',
                  'INPUT=' + bam_filename,
                  'OUTPUT=' + metrics_filename,
                  'REFERENCE_SEQUENCE=' + ref_genome,
                  'S=' + summary_filename,
                  'CHART_OUTPUT=' + chart_filename,
        'VALIDATION_STRINGENCY=LENIENT',
                  'TMP_DIR=' + tempdir,
        'MAX_RECORDS_IN_RAM=150000',
        docker_image=docker_image
    )


def bam_collect_insert_metrics(bam_filename, flagstat_metrics_filename,
                               metrics_filename, histogram_filename, tempdir,
                               mem="2G", docker_image=None):
    # Check if any paired reads exist
    has_paired = None
    with open(flagstat_metrics_filename) as f:
        for line in f:
            if 'properly paired' in line:
                if line.startswith('0 '):
                    has_paired = False
                else:
                    has_paired = True

    if has_paired is None:
        raise Exception('Unable to determine number of properly paired reads from {}'.format(
            flagstat_metrics_filename))

    if not has_paired:
        with open(metrics_filename, 'w') as f:
            f.write('## FAILED: No properly paired reads\n')
        with open(histogram_filename, 'w'):
            pass
        return

    if not os.path.exists(tempdir):
        makedirs(tempdir)

    pypeliner.commandline.execute(
        'picard', '-Xmx' + mem, '-Xms' + mem,
        '-XX:ParallelGCThreads=1',
        'CollectInsertSizeMetrics',
                  'INPUT=' + bam_filename,
                  'OUTPUT=' + metrics_filename,
                  'HISTOGRAM_FILE=' + histogram_filename,
        'ASSUME_SORTED=True',
        'VALIDATION_STRINGENCY=LENIENT',
                  'TMP_DIR=' + tempdir,
        'MAX_RECORDS_IN_RAM=150000',
        docker_image=docker_image
    )
