'''
Pipeline for alignment and QC of NextSeq sequencing data.
@author: Adi Steif

Steps:
- Convert BCL files to demultiplexed FASTQ format
- Produce FastQC report
- Align FASTQ files with BWA (include read group tags)
- Sort BAM files with Picard
- Mark duplicate reads with Picard
- Index the final BAM files with Samtools
- Extract metrics with Samtools flagstat and Picard
- Generate metric table
- Output basic QC plots

Notes: does NOT include indel realignment!

'''
from pipelines.io import make_directory, make_parent_directory
from ruffus import *
from ruffus.ruffus_utility import CHECKSUM_FILE_TIMESTAMPS

import argparse
import os
import yaml
import shutil
import warnings

#=======================================================================================================================
# Read Command Line Input
#=======================================================================================================================
parser = argparse.ArgumentParser()

parser.add_argument('config_file',
                    help='''Path to yaml config file.''')

parser.add_argument('--num_cpus', type=int, default=1,
                    help='''Number of cpus to use for the analysis. If set to -1 then as many cpus as samples will
                    be used. Default is 1.''')

parser.add_argument('--mode', choices=['local', 'cluster', 'printout'], default='printout',
                    help='''Mode to run the pipeline in. local will run the pipeline on the compute it was launched
                    from. cluster will submit the jobs to a cluster using SGE. printout shows which tasks will be
                    run. default is printout.''')

parser.add_argument('--install_dir', default=None,
                    help='''Path to local installation files to override system defaults.''')

args = parser.parse_args()

#=======================================================================================================================
# Read Input Files
#=======================================================================================================================
try:
    with open(args.config_file) as file:
        config = yaml.load(file)
        
except IOError as e:
    print 'Unable to open config file: {0}'.format(args.config_file)

sample_sheet_file = os.path.join(config['nextseq_dir'], 'SampleSheet.csv')

try:
    with open(sample_sheet_file) as file:
        lines = [x.strip('\n').strip(',') for x in file.readlines()]
    
    run_id = [s.split(',')[1] for s in lines if 'Experiment Name,' in s][0]
    
    library_id = [s.split(',')[1] for s in lines if 'Description,' in s][0]
    
except IOError as e:
    print 'Unable to open file \'SampleSheet.csv\' in directory: {0}'.format(config['nextseq_dir'])

#=======================================================================================================================
# Scripts
#=======================================================================================================================
cwd = os.path.dirname(os.path.realpath(__file__))

bin_dir = os.path.join(cwd, 'bin')

run_fastqc_script = os.path.join(bin_dir, 'run_fastqc.sh')

run_bwa_script = os.path.join(bin_dir, 'run_bwa_paired_end.sh')

run_samtools_flagstat_script = os.path.join(bin_dir, 'run_samtools_flagstat.sh')

extract_metrics_table_script = os.path.join(bin_dir, 'extract_metrics_table.py')

plot_metrics_script = os.path.join(bin_dir, 'plot_metrics.py')

#=======================================================================================================================
# Set System Paths
#=======================================================================================================================
if args.install_dir is not None:
    os.environ['PATH'] = os.path.join(args.install_dir, 'bin') + ':' + os.environ['PATH']
    
    if 'LD_LIBRARY_PATH' in os.environ:
        os.environ['LD_LIBRARY_PATH'] = os.path.join(args.install_dir, 'lib') + ':' + os.environ['LD_LIBRARY_PATH']
    
    else:
        os.environ['LD_LIBRARY_PATH'] = os.path.join(args.install_dir, 'lib')

#=======================================================================================================================
# Pipeline
#=======================================================================================================================
def generate_fastq_file_pairs():
    '''
    Note: could modify this to read the 'Project' for each cell, and make the output fastq
    directory take this into account (e.g. fastq/DLP/sample*.fastq.gz). This will allow 
    users to put libraries for multiple projects on a single run.
    '''
    start_index = lines.index('[Data]')+2
    
    num_samples = len(lines[start_index:])
    
    for i, line in zip(range(num_samples), lines[start_index:]):
        sample_id = line.split(',')[0]
        
        fastq_file_1 = os.path.join(config['out_dir'], 'fastq', '{0}_S{1}_R1_001.fastq.gz'.format(sample_id, str(i+1)))
        fastq_file_2 = os.path.join(config['out_dir'], 'fastq', '{0}_S{1}_R2_001.fastq.gz'.format(sample_id, str(i+1)))
        
        yield [fastq_file_1, fastq_file_2]

def generate_fastq_file_list():
    fastq_list = []
    
    for files in generate_fastq_file_pairs():
        fastq_list.append(files)
    
    yield [config['nextseq_dir'], fastq_list]

def generate_fastq_files_qc():
    for files in generate_fastq_file_pairs():
        fastq_file_1, fastq_file_2 = files
        
        sample_id = os.path.basename(fastq_file_1).split('_')[0]
        
        fastqc_html_1 = os.path.join(config['out_dir'], 'metrics', 'fastqc', '{0}_R1.fastqc.html'.format(sample_id))
            
        fastqc_zip_1 = os.path.join(config['out_dir'], 'metrics', 'fastqc', '{0}_R1.fastqc.zip'.format(sample_id))
            
        yield [[fastq_file_1], [fastqc_html_1, fastqc_zip_1]]
        
        fastqc_html_2 = os.path.join(config['out_dir'], 'metrics', 'fastqc', '{0}_R2.fastqc.html'.format(sample_id))
        
        fastqc_zip_2 = os.path.join(config['out_dir'], 'metrics', 'fastqc', '{0}_R2.fastqc.zip'.format(sample_id))

        yield [[fastq_file_2], [fastqc_html_2, fastqc_zip_2]]

def generate_fastq_files_align():
    for files in generate_fastq_file_pairs():
        fastq_file_1, fastq_file_2 = files
        
        sample_id = os.path.basename(fastq_file_1).split('_')[0]
        
        bam_file = os.path.join(config['out_dir'], 'tmp', '{0}.bam'.format(sample_id))

        yield [[fastq_file_1, fastq_file_2], bam_file]

@files(generate_fastq_file_list)
def demultiplex_fastq_files(nextseq_dir, out_files):
    out_dir = os.path.join(config['out_dir'], 'fastq')
    
    make_parent_directory(out_dir)
    
    cmd = 'bcl2fastq'
    
    cmd_args = [
                '--runfolder-dir=' + nextseq_dir, 
                '--output-dir=' + out_dir, 
                '--no-lane-splitting'
                ]
    
    run_cmd(cmd, cmd_args, max_mem=34)

@follows(demultiplex_fastq_files)
@files(generate_fastq_files_qc)
def produce_fastqc_report(in_file, out_files):
    make_parent_directory(out_files[0])
    
    cmd = 'sh'
    
    cmd_args = [
                run_fastqc_script, 
                in_file[0], 
                out_files[0], 
                out_files[1]
                ]
    
    run_cmd(cmd, cmd_args)

@follows(produce_fastqc_report)
@files(generate_fastq_files_align)
def align_fastq_files(in_files, out_file):
    make_parent_directory(out_file)
    
    tmp_file = out_file + '.tmp'
    
    cmd = 'sh'
    
    cmd_args = [
                run_bwa_script,
                config['ref_genome'],
                in_files[0],
                in_files[1],
                tmp_file
                ]
    
    if 'read_group' in config.keys():
        sample_id = os.path.basename(out_file).split('.')[0]
        
        read_group = ('@RG\tID:' + str(library_id) + '_' 
                                 + str(sample_id) + '_' 
                                 + str(run_id) + 
                         '\tPL:' + str(config['read_group']['PL']) +
                         '\tPU:' + str(run_id) +
                         '\tLB:' + str(library_id) + '_' + str(sample_id) +
                         '\tSM:' + str(sample_id) +
                         '\tCN:' + str(config['read_group']['CN']))
        
        cmd_args.append(read_group)
    
    else:
        warnings.warn('Config file does not contain read group information! ' + 
                      'This will affect duplicate marking if BAMs are later merged. ' +
                      'Creating BAM without read group information in header.')
    
    run_cmd(cmd, cmd_args, mem=10, max_mem=20)
    
    shutil.move(tmp_file, out_file)

@transform(align_fastq_files, suffix('.bam'), '.sorted.bam')
def sort_bam_file(bam_file, sorted_bam_file):
    tmp_file = sorted_bam_file + '.tmp'
    
    cmd = 'java'
    
    cmd_args = [
                '-Xmx16g',
                '-jar',
                config['picard_jar'],
                'SortSam',
                'INPUT=' + bam_file, 
                'OUTPUT=' + tmp_file,
                'SORT_ORDER=coordinate',
                'VALIDATION_STRINGENCY=LENIENT',
                'MAX_RECORDS_IN_RAM=5000000'
                ]
    
    run_cmd(cmd, cmd_args, max_mem=20)
    
    shutil.move(tmp_file, sorted_bam_file)

@transform(sort_bam_file, regex('(.*)/tmp/(.*)\.sorted.bam'), r'\1/bam/\2.sorted.markdups.bam')
def markdups_bam_file(bam_file, markdups_bam_file):
    make_parent_directory(markdups_bam_file)
    
    sample_id = os.path.basename(bam_file).split('.')[0]
    
    metrics_file = os.path.join(config['out_dir'], 
                                'metrics', 
                                'duplication_metrics', 
                                '{0}.markdups.txt'.format(sample_id))
    
    make_parent_directory(metrics_file)
    
    tmp_file = markdups_bam_file + '.tmp'
    
    cmd = 'java'
    
    cmd_args = [
                '-Xmx4g',
                '-jar',
                config['picard_jar'],
                'MarkDuplicates', 
                'INPUT=' + bam_file, 
                'OUTPUT=' + tmp_file, 
                'METRICS_FILE=' + metrics_file,
                'REMOVE_DUPLICATES=False', 
                'ASSUME_SORTED=True',
                'VALIDATION_STRINGENCY=LENIENT'
                ]
    
    run_cmd(cmd, cmd_args, max_mem=20)
    
    shutil.move(tmp_file, markdups_bam_file)

@transform(markdups_bam_file, suffix('.bam'), '.bam.bai')
def index_bam_file(bam_file, bai_file):
    cmd = 'samtools'
    
    cmd_args = ['index', bam_file]
    
    run_cmd(cmd, cmd_args)

@follows(index_bam_file)
@transform(markdups_bam_file, regex(r'(.*)/bam/(.*)\.sorted\.markdups\.bam'), 
                                    r'\1/metrics/flagstat_metrics/\2.flagstat_metrics.txt')
def extract_flagstat_metrics(bam_file, flagstat_file):
    make_parent_directory(flagstat_file)
    
    cmd = 'sh'
    
    cmd_args = [
                run_samtools_flagstat_script, 
                bam_file, 
                flagstat_file
                ]
    
    run_cmd(cmd, cmd_args)


@follows(index_bam_file)
@transform(markdups_bam_file, regex(r'(.*)/bam/(.*)\.sorted\.markdups\.bam'), 
                                    r'\1/metrics/wgs_metrics/\2.wgs_metrics.txt')
def extract_wgs_metrics(bam_file, metrics_file):
    make_parent_directory(metrics_file)
    
    # note this Picard command doesn't run if the ref genome .fai file
    # is in the same directory as the reference genome...
    #ref_genome = config['ref_genome'].replace('ref_genomes', 'ref_picard')
    #jar_file = os.path.join(config['picard_dir'], 'CollectWgsMetrics.jar')
    
    cmd = 'java'
    
    cmd_args = [
                '-Xmx4g',
                '-jar',
                config['picard_jar'],
                'CollectWgsMetrics', 
                'INPUT=' + bam_file, 
                'OUTPUT=' + metrics_file, 
                'REFERENCE_SEQUENCE=' + config['ref_genome'], 
                'MINIMUM_BASE_QUALITY=' + str(config['min_bqual']), 
                'MINIMUM_MAPPING_QUALITY=' + str(config['min_mqual']), 
                'COVERAGE_CAP=500', 
                'VALIDATION_STRINGENCY=LENIENT'
                ]
    
    if 'count_unpaired' in config.keys():
        cmd_args.append('COUNT_UNPAIRED=True')
    
    run_cmd(cmd, cmd_args, max_mem=10)

@follows(index_bam_file)
@transform(markdups_bam_file, regex(r'(.*)/bam/(.*)\.sorted\.markdups\.bam'), 
                                    r'\1/metrics/insert_metrics/\2.insert_metrics.txt')
def extract_insert_metrics(bam_file, metrics_file):
    make_parent_directory(metrics_file)
    
    hist_file = metrics_file.replace('.txt', '.pdf')
    
    cmd = 'java'
    
    cmd_args = [
                '-Xmx4g',
                '-jar',
                config['picard_jar'], 
                'CollectInsertSizeMetrics',
                'INPUT=' + bam_file, 
                'OUTPUT=' + metrics_file, 
                'HISTOGRAM_FILE=' + hist_file, 
                'ASSUME_SORTED=True',
                'VALIDATION_STRINGENCY=LENIENT'
                ]
    
    run_cmd(cmd, cmd_args, max_mem=10)

@follows(extract_flagstat_metrics, extract_wgs_metrics, extract_insert_metrics)
@collate(markdups_bam_file, regex(r'(.*)/bam/.*bam'), 
                                  r'\1/metrics/summary/{0}_{1}.metrics.csv'.format(library_id, run_id))
def extract_metrics_table(bam_files, out_file):
    make_parent_directory(out_file)
    
    cmd = 'python'
    
    cmd_args = [
                extract_metrics_table_script,
                os.path.join(config['out_dir'], 'metrics'),
                out_file
                ]
    
    run_cmd(cmd, cmd_args)

@transform(extract_metrics_table, suffix('.csv'), '.pdf')
def plot_metrics(metrics_table, out_file):
    make_parent_directory(out_file)
    
    cmd = 'python'
    
    cmd_args = [
                plot_metrics_script,
                sample_sheet_file, 
                metrics_table,
                out_file
                ]
    
    run_cmd(cmd, cmd_args)

@follows(plot_metrics)
def end():
    pass

#=======================================================================================================================
# Run Pipeline
#=======================================================================================================================    
if args.mode in ['cluster', 'local']:
    if args.mode == 'cluster':
        from pipelines.job_manager import ClusterJobManager
        
        import datetime
        
        log_dir = os.path.join(config['out_dir'], 'log', datetime.datetime.now().isoformat(','))
        
        job_manager = ClusterJobManager(log_dir)
    
    elif args.mode == 'local':
        from pipelines.job_manager import LocalJobManager
        
        job_manager = LocalJobManager()
    
    run_cmd = job_manager.run_job
    
    try:
        pipeline_run(end, multithread=args.num_cpus, checksum_level=CHECKSUM_FILE_TIMESTAMPS)
    
    finally:
        job_manager.close()

elif args.mode == 'printout':
    import sys
    
    pipeline_printout(sys.stdout, end, verbose=3, wrap_width=200)
