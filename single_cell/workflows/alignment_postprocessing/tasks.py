'''
Created on Jul 24, 2017

@author: dgrewal
'''
import pypeliner
from scripts import GenerateCNMatrix
from scripts import CollectMetrics
import pandas as pd
import os
import errno



def makedirs(directory):
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def postprocess_bam(infile, outfile, outfile_index, tempdir,
                    config, markdups_metrics, flagstat_metrics):
    
    if not os.path.exists(tempdir):
        makedirs(tempdir)
    
    sorted_bam = os.path.join(tempdir, 'sorted.bam')
    
    bam_sort(infile, sorted_bam, tempdir)

    bam_markdups(sorted_bam, outfile, markdups_metrics, tempdir)

    bam_index(outfile, outfile_index)
    
    bam_flagstat(outfile, flagstat_metrics)
    

def bam_flagstat(infile, outfile):
    pypeliner.commandline.execute(
        'samtools', 'flagstat',
        infile, '>',
        outfile
        )


def bam_index(infile, outfile):
    pypeliner.commandline.execute(
        'samtools', 'index',
        infile,
        outfile
        )

def bam_sort(bam_filename, sorted_bam_filename, tempdir):
    if not os.path.exists(tempdir):
        makedirs(tempdir)

    pypeliner.commandline.execute(
        'picard', '-Xmx2G', '-Xms2G',
        '-XX:ParallelGCThreads=1',
        'SortSam',
        'INPUT=' + bam_filename,
        'OUTPUT=' + sorted_bam_filename,
        'SORT_ORDER=coordinate',
        'VALIDATION_STRINGENCY=LENIENT',
        'TMP_DIR='+tempdir,
        'MAX_RECORDS_IN_RAM=150000')


def bam_markdups(bam_filename, markduped_bam_filename, metrics_filename, tempdir):
    if not os.path.exists(tempdir):
        makedirs(tempdir)

    pypeliner.commandline.execute(
        'picard', '-Xmx2G', '-Xms2G',
        '-XX:ParallelGCThreads=1',
        'MarkDuplicates',
        'INPUT=' + bam_filename,
        'OUTPUT=' + markduped_bam_filename,
        'METRICS_FILE=' + metrics_filename,
        'REMOVE_DUPLICATES=False',
        'ASSUME_SORTED=True',
        'VALIDATION_STRINGENCY=LENIENT',
        'TMP_DIR='+tempdir,
        'MAX_RECORDS_IN_RAM=150000')


def bam_collect_wgs_metrics(bam_filename, ref_genome, metrics_filename, config, tempdir):

    if not os.path.exists(tempdir):
        makedirs(tempdir)

    pypeliner.commandline.execute(
        'picard', '-Xmx2G', '-Xms2G',
        '-XX:ParallelGCThreads=1',
        'CollectWgsMetrics',
        'INPUT=' + bam_filename,
        'OUTPUT=' + metrics_filename,
        'REFERENCE_SEQUENCE=' + ref_genome,
        'MINIMUM_BASE_QUALITY=' + str(config['picard_wgs_params']['min_bqual']),
        'MINIMUM_MAPPING_QUALITY=' + str(config['picard_wgs_params']['min_mqual']),
        'COVERAGE_CAP=500',
        'VALIDATION_STRINGENCY=LENIENT',
        'COUNT_UNPAIRED=' + ('True' if config['picard_wgs_params']['count_unpaired'] else 'False'),
        'TMP_DIR='+tempdir,
        'MAX_RECORDS_IN_RAM=150000')


def bam_collect_gc_metrics(bam_filename, ref_genome, metrics_filename, summary_filename, chart_filename, tempdir):
    if not os.path.exists(tempdir):
        makedirs(tempdir)

    pypeliner.commandline.execute(
        'picard', '-Xmx2G', '-Xms2G',
        '-XX:ParallelGCThreads=1',
        'CollectGcBiasMetrics',
        'INPUT=' + bam_filename,
        'OUTPUT=' + metrics_filename,
        'REFERENCE_SEQUENCE=' + ref_genome,
        'S=' + summary_filename,
        'CHART_OUTPUT=' + chart_filename,
        'VALIDATION_STRINGENCY=LENIENT',
        'TMP_DIR='+tempdir,
        'MAX_RECORDS_IN_RAM=150000')

def bam_collect_insert_metrics(bam_filename, flagstat_metrics_filename, metrics_filename, histogram_filename, tempdir):
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
        raise Exception('Unable to determine number of properly paired reads from {}'.format(flagstat_metrics_filename))

    if not has_paired:
        with open(metrics_filename, 'w') as f:
            f.write('## FAILED: No properly paired reads\n')
        with open(histogram_filename, 'w'):
            pass
        return

    if not os.path.exists(tempdir):
        makedirs(tempdir)

    pypeliner.commandline.execute(
        'picard', '-Xmx2G', '-Xms2G',
        '-XX:ParallelGCThreads=1',
        'CollectInsertSizeMetrics',
        'INPUT=' + bam_filename,
        'OUTPUT=' + metrics_filename,
        'HISTOGRAM_FILE=' + histogram_filename,
        'ASSUME_SORTED=True',
        'VALIDATION_STRINGENCY=LENIENT',
        'TMP_DIR='+tempdir,
        'MAX_RECORDS_IN_RAM=150000')

def collect_gc(infile, outfile, sample_id):
    gen_gc = GenerateCNMatrix(infile, outfile, ',',
                              'NORMALIZED_COVERAGE', sample_id,
                              'gcbias')
    gen_gc.main()

def collect_metrics(flagstat_metrics, markdups_metrics, insert_metrics,
                    wgs_metrics, output, merged_metrics):

    for sample in flagstat_metrics.keys():
        flgstat = flagstat_metrics[sample]
        mkdup = markdups_metrics[sample]
        insrt = insert_metrics[sample]
        wgs = wgs_metrics[sample]
        outfile = output[sample]

        collmet = CollectMetrics(wgs, insrt, flgstat,
                                 mkdup, outfile, sample)
        collmet.main()

    samples = flagstat_metrics.keys()
    concatenate_csv(output, merged_metrics, samples)

def concatenate_csv(in_filenames, out_filename, sample_ids, nan_val='NA'):
    data = []

    if not isinstance(in_filenames, dict):
        in_filenames = dict([(sampid, in_filenames[sampid])
                             for sampid in sample_ids])
    
    for _, in_filename in in_filenames.iteritems():
        with open(in_filename) as f:
            first_line = f.readline()
            if len(first_line) == 0:
                continue
        data.append(pd.read_csv(in_filename, dtype=str))
    data = pd.concat(data, ignore_index=True)
    data = data.fillna(nan_val)
    data.to_csv(out_filename, index=False)


def merge_csv(in_filenames, out_filename, how, on, sample_ids, nan_val='NA'):
    data = []

    if not isinstance(in_filenames, dict):
        in_filenames = dict([(sampid, in_filenames[sampid])
                             for sampid in sample_ids])

    
    for _, in_filename in in_filenames.iteritems():
        with open(in_filename) as f:
            first_line = f.readline()
            if len(first_line) == 0:
                continue
        data.append(pd.read_csv(in_filename, dtype=str))

    data = merge_frames(data, how, on)
    data = data.fillna(nan_val)
    data.to_csv(out_filename, index=False)


def merge_frames(frames, how, on):
    '''
    annotates input_df using ref_df
    '''

    if ',' in on:
        on = on.split(',')

    if len(frames) == 1:
        return frames[0]
    else:
        left = frames[0]
        right = frames[1]
        merged_frame = pd.merge(left, right,
                                how=how,
                                on=on)
        for frame in frames[2:]:
            merged_frame = pd.merge(merged_frame, frame,
                                    how=how,
                                    on=on)
        return merged_frame
