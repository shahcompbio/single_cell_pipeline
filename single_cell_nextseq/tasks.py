import os
import errno
import pandas as pd
import pypeliner.commandline


def makedirs(directory):
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def generate_fastq_file_pairs(sample_sheet_filename, fastq_directory):
    ''' Generate standardized filenames from sample sheet.
    '''
    with open(sample_sheet_filename) as f:
        lines = [x.strip('\n').strip(',') for x in f.readlines()]

    start_index = lines.index('[Data]')+2

    num_samples = len(lines[start_index:])

    for i, line in zip(range(num_samples), lines[start_index:]):
        sample_id = line.split(',')[0]

        fastq_file_1 = os.path.join(fastq_directory, '{0}_S{1}_R1_001.fastq.gz'.format(sample_id, str(i+1)))
        fastq_file_2 = os.path.join(fastq_directory, '{0}_S{1}_R2_001.fastq.gz'.format(sample_id, str(i+1)))

        yield sample_id, fastq_file_1, fastq_file_2


def demultiplex_fastq_files(sample_sheet_filename, nextseq_directory, fastq_1_filenames, fastq_2_filenames, temp_directory):
    pypeliner.commandline.execute(
        'bcl2fastq',
        '--runfolder-dir=' + nextseq_directory,
        '--output-dir=' + temp_directory,
        '--no-lane-splitting')

    for sample_id, temp_fq_1, temp_fq_2 in generate_fastq_file_pairs(sample_sheet_filename, temp_directory):
        os.rename(temp_fq_1, fastq_1_filenames[sample_id])
        os.rename(temp_fq_2, fastq_2_filenames[sample_id])


def parse_sample_sheet(sample_sheet_filename, sample_info_filename):
    """ Parse the sample sheet data into a csv table
    """
    header = None
    with open(sample_sheet_filename) as sample_sheet_file:
        for idx, line in enumerate(sample_sheet_file):
            if line.startswith('[Data]'):
                header = idx + 1
                break

    if header is None:
        raise Exception('Unable to find [Data] in samplesheet {}'.format(sample_sheet_filename))

    data = pd.read_csv(sample_sheet_filename, header=header, dtype=str)

    column_renames = {
        'Sample_ID': 'sample_id',
        'Sample_Well': 'sample_well',
        'Description': 'sample_description',
        'Sample_Plate': 'sample_plate',
        'I5_Index_ID': 'i5_barcode',
        'I7_Index_ID': 'i7_barcode'}
    
    data = data.rename(columns=column_renames)
    data = data[column_renames.values()]

    def extract_cell_call(description):
        if ';' in description:
            return description.split(';')[0].replace('CC=','')
        elif description == '':
            return 'C1'
        else:
            return description

    data['cell_call'] = data['sample_description'].apply(extract_cell_call)

    def extract_experimental_condition(description):
        if ';' in description:
            return description.split(';')[1].replace('EC=','')
        elif description == '':
            return 'A'
        else:
            return 'NA'

    data['experimental_condition'] = data['sample_description'].apply(extract_experimental_condition)

    # Set to default values
    data.loc[data['sample_well'] == '', 'sample_well'] = 'R1_C1'
    data.loc[data['sample_plate'] == '', 'sample_plate'] = 'R1-C1'
    data.loc[data['i5_barcode'] == '', 'i5_barcode'] = 'i5-1'
    data.loc[data['i7_barcode'] == '', 'i7_barcode'] = 'i7-1'

    # Replace '-' with '_' for sample id and sample plate
    data['sample_id'] = data['sample_id'].apply(lambda a: a.replace('-', '_'))
    data['sample_plate'] = data['sample_plate'].apply(lambda a: a.replace('-', '_'))

    data.to_csv(sample_info_filename, index=False)


def produce_fastqc_report(fastq_filename, fastq_basename, output_html, output_plots, temp_directory):
    makedirs(temp_directory)

    pypeliner.commandline.execute(
        'fastqc',
        '--outdir=' + temp_directory,
        fastq_filename)

    output_basename = os.path.join(temp_directory, fastq_basename)
    os.rename(output_basename + '_fastqc.html', output_html)
    os.rename(output_basename + '_fastqc.zip', output_plots)


def run_trimgalore(fastq_1_filename, fastq_2_filename, fastq_1_basename, fastq_2_basename,
                   trim_1_filename, trim_2_filename, results_directory, adapter, adapter2):
    makedirs(results_directory)

    print ' '.join((
        'trim_galore',
        '--fastqc',
        '--paired',
        '--output_dir', results_directory,
        '--adapter', adapter,
        '--adapter2', adapter2,
        fastq_1_filename,
        fastq_2_filename))

    pypeliner.commandline.execute(
        'trim_galore',
        '--fastqc',
        '--paired',
        '--output_dir', results_directory,
        '--adapter', adapter,
        '--adapter2', adapter2,
        fastq_1_filename,
        fastq_2_filename)

    os.rename(os.path.join(results_directory, fastq_1_basename) + '_val_1.fq.gz', trim_1_filename)
    os.rename(os.path.join(results_directory, fastq_2_basename) + '_val_2.fq.gz', trim_2_filename)


def bam_sort(bam_filename, sorted_bam_filename):
    pypeliner.commandline.execute(
        'picard', '-Xmx12G',
        'SortSam',
        'INPUT=' + bam_filename,
        'OUTPUT=' + sorted_bam_filename,
        'SORT_ORDER=coordinate',
        'VALIDATION_STRINGENCY=LENIENT',
        'MAX_RECORDS_IN_RAM=5000000')


def bam_markdups(bam_filename, markduped_bam_filename, metrics_filename):
    pypeliner.commandline.execute(
        'picard', '-Xmx8G',
        'MarkDuplicates',
        'INPUT=' + bam_filename,
        'OUTPUT=' + markduped_bam_filename,
        'METRICS_FILE=' + metrics_filename,
        'REMOVE_DUPLICATES=False',
        'ASSUME_SORTED=True',
        'VALIDATION_STRINGENCY=LENIENT')


def bam_collect_wgs_metrics(bam_filename, ref_genome, metrics_filename, config):
    pypeliner.commandline.execute(
        'picard', '-Xmx12G',
        'CollectWgsMetrics',
        'INPUT=' + bam_filename,
        'OUTPUT=' + metrics_filename,
        'REFERENCE_SEQUENCE=' + ref_genome,
        'MINIMUM_BASE_QUALITY=' + str(config['min_bqual']),
        'MINIMUM_MAPPING_QUALITY=' + str(config['min_mqual']),
        'COVERAGE_CAP=500',
        'VALIDATION_STRINGENCY=LENIENT',
        'COUNT_UNPAIRED=' + ('True' if config['count_unpaired'] else 'False'))


def bam_collect_gc_metrics(bam_filename, ref_genome, metrics_filename, summary_filename, chart_filename):
    pypeliner.commandline.execute(
        'picard', '-Xmx12G',
        'CollectGcBiasMetrics',
        'INPUT=' + bam_filename,
        'OUTPUT=' + metrics_filename,
        'REFERENCE_SEQUENCE=' + ref_genome,
        'S=' + summary_filename,
        'CHART_OUTPUT=' + chart_filename,
        'VALIDATION_STRINGENCY=LENIENT')


def bam_collect_insert_metrics(bam_filename, flagstat_metrics_filename, metrics_filename, histogram_filename):
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

    pypeliner.commandline.execute(
        'picard', '-Xmx12G',
        'CollectInsertSizeMetrics',
        'INPUT=' + bam_filename,
        'OUTPUT=' + metrics_filename,
        'HISTOGRAM_FILE=' + histogram_filename,
        'ASSUME_SORTED=True',
        'VALIDATION_STRINGENCY=LENIENT')


scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
run_hmmcopy_rscript = os.path.join(scripts_directory, 'hmmcopy.R')


def run_hmmcopy(
    readcount_wig_filename,
    corrected_reads_filename,
    segments_filename,
    parameters_filename,
    posterior_marginals_filename,
    sample_id,
    config):

    pypeliner.commandline.execute(
        'Rscript', run_hmmcopy_rscript,
        '--tumour_file=' + readcount_wig_filename,
        '--gc_file=' + config['gc_wig_file'],
        '--map_file=' + config['map_wig_file'],
        '--reads_output=' + corrected_reads_filename,
        '--segs_output=' + segments_filename,
        '--params_output=' + parameters_filename,
        '--post_marginals_output=' + posterior_marginals_filename,
        '--map_cutoff=' + str(config['map_cutoff']),
        '--num_states=' + str(config['num_states']),
        '--param_mu=' + str(config['parameters']['mu']),
        '--param_m=' + str(config['parameters']['m']),
        '--param_k=' + str(config['parameters']['kappa']),
        '--param_e=' + str(config['parameters']['e']),
        '--param_s=' + str(config['parameters']['s']),
        '--sample_id=' + sample_id)


def concatenate_csv(in_filenames, out_filename):
    data = []
    for key, in_filename in in_filenames.iteritems():
        with open(in_filename) as f:
            first_line = f.readline()
            if len(first_line) == 0:
                continue
        data.append(pd.read_csv(in_filename, dtype=str))
    data = pd.concat(data, ignore_index=True)
    data.to_csv(out_filename, index=False)

