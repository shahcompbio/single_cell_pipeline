'''
Pipeline for alignment and QC of NextSeq sequencing data.
@author: Adi Steif

Inputs:
- Config file in .yaml format
    - Includes path to NextSeq directory, which must contain a SampleSheet.csv file
    - Note SampleSheet.csv must follow specific conventions beyond Illumina's requirements

Steps:
- Convert BCl files to demultiplexed FASTQ format
- Generate FastQC report
- Align FASTQ files with BWA (include read group tags)
- Sort BAM files with Picard
- Mark duplicate reads with Picard
- Index the final BAM files with Samtools
- Extract metrics with Samtools flagstat and Picard
- Generate metric table
- Output basic QC plots

Requirements:
- Python 2.7
    - Andy Roth's 'pipelines' package
    - ruffus
    - drmaa
    - yaml
    - pandas
    - matplotlib
    - seaborn
- Illumina's bcl2fastq
- FastQC
- Samtools
- BWA
- Picard
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
    start_index = lines.index('[Data]')+2
    
    num_samples = len(lines[start_index:])
    
    fastq_dir = os.path.join(config['nextseq_dir'], 'Data', 'Intensities', 'BaseCalls')
    
    for i, line in zip(range(num_samples), lines[start_index:]):
        sample_id = line.split(',')[0]
        
        fastq_file_1 = os.path.join(fastq_dir, '{0}_S{1}_R1_001.fastq.gz'.format(sample_id, str(i+1)))
        fastq_file_2 = os.path.join(fastq_dir, '{0}_S{1}_R2_001.fastq.gz'.format(sample_id, str(i+1)))
        
        yield [fastq_file_1, fastq_file_2]

@originate(generate_fastq_file_pairs)
def demultiplex_fastq_files(out_files):
    cmd = 'bcl2fastq'
    
    cmd_args = [
                '--no-lane-splitting', 
                '--runfolder-dir ' + config['nextseq_dir'], 
                ]
    
    run_cmd(cmd, cmd_args)

@subdivide(demultiplex_fastq_files, regex(r'.*/BaseCalls/(.*)_S[0-9]+\_R(\d)\_001\.fastq\.gz'), 
                                         [r'{0}/metrics/fastqc/\1_R\2.fastqc.html'.format(config['out_dir']), 
                                          r'{0}/metrics/fastqc/\1_R\2.fastqc.zip'.format(config['out_dir'])])
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

@transform(demultiplex_fastq_files, regex(r'.*/BaseCalls/(.*)_S[0-9]+\_R\d\_001\.fastq\.gz'), 
                                          r'{0}/tmp/\1.bam'.format(config['out_dir']))
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
    jar_file = os.path.join(config['picard_dir'], 'picard.jar')
    
    tmp_file = sorted_bam_file + '.tmp'
    
    cmd = 'java'
    
    cmd_args = [
                '-Xmx16g',
                '-jar',
                jar_file,
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
    
    metrics_file = os.path.join(config['out_dir'], 'metrics', 'duplication_metrics', '{0}.markdups.txt'.format(sample_id))
    
    make_parent_directory(metrics_file)
    
    tmp_file = markdups_bam_file + '.tmp'
    
    jar_file = os.path.join(config['picard_dir'], 'picard.jar')
    
    cmd = 'java'
    
    cmd_args = [
                '-Xmx4g',
                '-jar',
                jar_file,
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
@collate([markdups_bam_file, rmdups_bam_file], regex('(.*)/bam/.*bam$'), r'\1/metrics/summary/{0}.metrics_table.csv'.format(config['analysis_id']))
def run_metrics_pipeline(bam_files, metrics_file):
    cmd = 'sh'
    
    cmd_args = [
                run_metrics_pipeline_script,
                args.install_dir,
                args.config_file, 
                args.num_cpus, 
                args.mode
                ]
    
    run_cmd(cmd, cmd_args)

@follows(run_metrics_pipeline)
def end():
    pass

#=======================================================================================================================
# Run pipeline
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
        # pipeline_run(end, multiprocess=args.num_cpus, use_multi_threading=True)
        pipeline_run(end, multithread=args.num_cpus, checksum_level=CHECKSUM_FILE_TIMESTAMPS)
    
    finally:
        job_manager.close()

elif args.mode == 'printout':
    import sys
    
    pipeline_printout(sys.stdout, end, verbose=3, wrap_width=200)
