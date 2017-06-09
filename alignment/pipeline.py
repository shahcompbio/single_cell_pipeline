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
from glob import glob

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

run_read_counter_script = os.path.join(bin_dir, 'run_read_counter.sh')

run_insert_metrics_script = os.path.join(bin_dir, 'run_insert_metrics.sh')
run_gc_metrics_script = os.path.join(bin_dir, 'run_gc_metrics.sh')

run_hmmcopy_rscript = os.path.join(bin_dir, 'run_hmmcopy_single_cell.R')

run_hmmcopy_shscript = os.path.join(bin_dir, 'run_hmmcopy.sh')

extract_quality_metrics_script = os.path.join(bin_dir, 'extract_quality_metrics.py')

plot_hmmcopy_script = os.path.join(bin_dir, 'plot_hmmcopy.py')
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
def get_metrics_files():
    files = '{0}/metrics/summary/(.*)metrics.csv'.format(config['out_dir'])
    
    files = glob(files)

    out = '{0}/metrics/summary/allmetrics.csv'.format(config['out_dir'])

    print files,out

    yield [files,out]


def load_counts_files():
    tmp_dir = os.path.join(config['out_dir'], 'tmp')
    
    for file in os.listdir(config['bam_dir']):
        if file.endswith(".markdups.bam"):
            sample_id = file.split(".")[0]
            
            tumour_bam_file = os.path.join(config['bam_dir'], file)
            
            yield tumour_bam_file, os.path.join(tmp_dir, '{0}.tumour.wig'.format(sample_id))    

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
    
    cmd = config['bcl2fastq']
    
    cmd_args = [
                '--runfolder-dir=' + nextseq_dir, 
                '--output-dir=' + out_dir, 
                '--no-lane-splitting'
                ]
    
    run_cmd(cmd, cmd_args, max_mem=34, queue=config['queue'], hosts=config['hosts'])

@follows(demultiplex_fastq_files)
@files(generate_fastq_files_qc)
def produce_fastqc_report(in_file, out_files):
    make_parent_directory(out_files[0])
    
    cmd = 'sh'
    
    cmd_args = [
                run_fastqc_script, 
                in_file[0], 
                out_files[0], 
                out_files[1],
                config['fastqc'],
                config['java']
                ]
    
    run_cmd(cmd, cmd_args, queue=config['queue'], hosts=config['hosts'], max_mem=10)

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
                tmp_file,
                config['bwa'],
                config['samtools']
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

    run_cmd(cmd, cmd_args, mem=10, max_mem=20, queue=config['queue'], hosts=config['hosts'])
    
    shutil.move(tmp_file, out_file)

@transform(align_fastq_files, suffix('.bam'), '.sorted.bam')
def sort_bam_file(bam_file, sorted_bam_file):
    tmp_file = sorted_bam_file + '.tmp'
    
    cmd = config['java']
    
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
    
    run_cmd(cmd, cmd_args, max_mem=30, queue=config['queue'], hosts=config['hosts'])
    
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
    
    cmd = config['java']
    
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
    
    run_cmd(cmd, cmd_args, max_mem=20, queue=config['queue'], hosts=config['hosts'])
    
    shutil.move(tmp_file, markdups_bam_file)

@transform(markdups_bam_file, suffix('.bam'), '.bam.bai')
def index_bam_file(bam_file, bai_file):
    cmd = config['samtools']
    
    cmd_args = ['index', bam_file]
    
    run_cmd(cmd, cmd_args, queue=config['queue'], hosts=config['hosts'])

@follows(index_bam_file)
@transform(markdups_bam_file, regex(r'(.*)/bam/(.*)\.sorted\.markdups\.bam'), 
                                    r'\1/metrics/flagstat_metrics/\2.flagstat_metrics.txt')
def extract_flagstat_metrics(bam_file, flagstat_file):
    make_parent_directory(flagstat_file)
    
    cmd = 'sh'
    
    cmd_args = [
                run_samtools_flagstat_script, 
                bam_file, 
                flagstat_file,
                config['samtools']
                ]
    
    run_cmd(cmd, cmd_args, queue=config['queue'], hosts=config['hosts'])


@follows(index_bam_file)
@transform(markdups_bam_file, regex(r'(.*)/bam/(.*)\.sorted\.markdups\.bam'), 
                                    r'\1/metrics/wgs_metrics/\2.wgs_metrics.txt')
def extract_wgs_metrics(bam_file, metrics_file):
    make_parent_directory(metrics_file)
    
    # note this Picard command doesn't run if the ref genome .fai file
    # is in the same directory as the reference genome...
    #ref_genome = config['ref_genome'].replace('ref_genomes', 'ref_picard')
    #jar_file = os.path.join(config['picard_dir'], 'CollectWgsMetrics.jar')
    
    cmd = config['java']
    
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
    
    run_cmd(cmd, cmd_args, max_mem=10, queue=config['queue'], hosts=config['hosts'])


@follows(index_bam_file)
@transform(markdups_bam_file, regex(r'(.*)/bam/(.*)\.sorted\.markdups\.bam'), 
                                    r'\1/metrics/gc_metrics/\2.gc_metrics.txt')
def extract_gc_metrics(bam_file, metrics_file):
    make_parent_directory(metrics_file)
    
    # note this Picard command doesn't run if the ref genome .fai file
    # is in the same directory as the reference genome...
    #ref_genome = config['ref_genome'].replace('ref_genomes', 'ref_picard')
    #jar_file = os.path.join(config['picard_dir'], 'CollectWgsMetrics.jar')
    
    sum_file = metrics_file.replace('.txt', '.summ.txt')
    chart_file = metrics_file.replace('.txt', '.pdf')


    cmd = 'sh'

    cmd_args = [
                run_gc_metrics_script,
                config['java'],
                config['picard_jar'],
                bam_file,
                metrics_file,
                config['ref_genome'],
                sum_file,
                chart_file,
                config['R']
                ]




    # cmd = config['java']
    
    # cmd_args = [
    #             '-Xmx4g',
    #             '-jar',
    #             config['picard_jar'],
    #             'CollectGcBiasMetrics', 
    #             'INPUT=' + bam_file, 
    #             'OUTPUT=' + metrics_file, 
    #             'REFERENCE_SEQUENCE=' + config['ref_genome'], 
    #             'S=' + sum_file, 
    #             'CHART_OUTPUT=' + chart_file, 
    #             'VALIDATION_STRINGENCY=LENIENT'
    #             ]
    
    run_cmd(cmd, cmd_args, max_mem=10, queue=config['queue'], hosts=config['hosts'])


@follows(index_bam_file)
@transform(markdups_bam_file, regex(r'(.*)/bam/(.*)\.sorted\.markdups\.bam'), 
                                    r'\1/metrics/insert_metrics/\2.insert_metrics.txt')
def extract_insert_metrics(bam_file, metrics_file):
    make_parent_directory(metrics_file)
    
    hist_file = metrics_file.replace('.txt', '.pdf')


    cmd = 'sh'

    cmd_args = [
                run_insert_metrics_script,
                config['java'],
                config['picard_jar'],
                bam_file,
                metrics_file,
                hist_file,
                config['R']
                ]
    
    # cmd = config['java']
    
    # cmd_args = [
    #             '-Xmx4g',
    #             '-jar',
    #             config['picard_jar'], 
    #             'CollectInsertSizeMetrics',
    #             'INPUT=' + bam_file, 
    #             'OUTPUT=' + metrics_file, 
    #             'HISTOGRAM_FILE=' + hist_file, 
    #             'ASSUME_SORTED=True',
    #             'VALIDATION_STRINGENCY=LENIENT'
    #             ]
    
    run_cmd(cmd, cmd_args, max_mem=10, queue=config['queue'], hosts=config['hosts'])

@follows(extract_flagstat_metrics, extract_wgs_metrics, extract_insert_metrics)
@collate(markdups_bam_file, regex(r'(.*)/bam/.*bam'), 
                                  r'\1/metrics/summary/{0}_{1}.metrics.csv'.format(library_id, run_id))
def extract_metrics_table(bam_files, out_file):
    make_parent_directory(out_file)
    
    cmd = config['python']
    
    cmd_args = [
                extract_metrics_table_script,
                os.path.join(config['out_dir'], 'metrics'),
                out_file,
                '--samplesheet', config['library_info']
                ]
    
    run_cmd(cmd, cmd_args, queue=config['queue'], hosts=config['hosts'])


@follows(index_bam_file)
@transform(markdups_bam_file, regex(r'(.*)/bam/(.*)\.sorted\.markdups\.bam'), 
                                    r'\1/tmp/\2.tumour.wig')
def build_counts_files(in_file, out_file):
    make_parent_directory(out_file)
    
    cmd = 'sh'

    chroms = get_wig_chromosomes(config['gc_wig_file'])
    
    cmd_args = [run_read_counter_script,
                config['readcounter'],
                in_file, 
                out_file,
                config['bin_size'],
                config['min_mqual'],
                ','.join(chroms),
                ]
    run_cmd(cmd, cmd_args, queue=config['queue'], hosts=config['hosts'])





@transform(build_counts_files, regex('(.*)/tmp/(.*)\.(normal|tumour)\.wig'), (r'\1/results/\2.corrected_reads.csv', 
                                                                              r'\1/results/\2.segments.csv',
                                                                              r'\1/results/\2.parameters.csv',
                                                                              r'\1/results/\2.posterior_marginals.csv'))
def run_hmmcopy(in_file, out_files):
    make_parent_directory(out_files[0])
    
    out_dir = os.path.dirname(out_files[0])
    
    out_basename = os.path.basename(out_files[0]).split('.')[0]
    

    cmd = 'sh'

    cmd_args=[run_hmmcopy_shscript,
              config['Rscript'],
              config['R_libs'],
              run_hmmcopy_rscript,
              in_file,
              config['gc_wig_file'],
              config['map_wig_file'],
              out_dir,
              out_basename,
              config['map_cutoff'],
              config['num_states'],
              config['parameters']['mu'],
              config['parameters']['m'],
              config['parameters']['kappa'],
              config['parameters']['e'],
              config['parameters']['s'],

             ]

    # cmd = config['Rscript']
    
    # cmd_args = [run_hmmcopy_script,
    #             '--tumour_file=' + in_file, 
    #             '--gc_file=' + config['gc_wig_file'], 
    #             '--map_file=' + config['map_wig_file'], 
    #             '--out_dir=' + out_dir, 
    #             '--out_basename=' + out_basename, 
    #             ]
    
    # if 'map_cutoff' in config:
    #         cmd_args.append('--map_cutoff=' + str(config['map_cutoff']))
    
    # if 'num_states' in config:
    #         cmd_args.append('--num_states=' + str(config['num_states']))
    
    # if 'parameters' in config:
    #     if 'mu' in config['parameters']:
    #         cmd_args.append('--param_mu=' + config['parameters']['mu'])
        
    #     if 'm' in config['parameters']:
    #         cmd_args.append('--param_m=' + config['parameters']['m'])
        
    #     if 'kappa' in config['parameters']:
    #         cmd_args.append('--param_k=' + config['parameters']['kappa'])
        
    #     if 'e' in config['parameters']:
    #         cmd_args.append('--param_e=' + str(config['parameters']['e']))
        
    #     if 'gamma' in config['parameters']:
    #         cmd_args.append('--param_g=' + str(config['parameters']['gamma']))
        
    #     if 's' in config['parameters']:
    #         cmd_args.append('--param_s=' + str(config['parameters']['s']))
    
    run_cmd(cmd, cmd_args, queue=config['queue'], hosts=config['hosts'])

@merge(run_hmmcopy, '{0}/metrics/{1}.quality_metrics.csv'.format(config['out_dir'], config['analysis_id']))
def extract_quality_metrics(in_files, out_file):
    make_parent_directory(out_file)
    
    results_dir = os.path.dirname(in_files[0][0])
    
    cmd = config['python']
    
    cmd_args = [
                extract_quality_metrics_script, 
                '--hmmcopy_results_dir', results_dir,
                '--out_file', out_file
                ]
    
    run_cmd(cmd, cmd_args, max_mem=20, queue=config['queue'], hosts=config['hosts'])



@merge(run_hmmcopy, '{0}/metrics/{1}.segments.csv'.format(config['out_dir'], config['analysis_id']))
def merge_segs(in_files, out_file):
    make_parent_directory(out_file)
    
    merge_segs_script = os.path.join(bin_dir, 'merge.py')

    # in_segfiles = glob(os.path.dirname(in_files[0][0]) + '/*.segments.csv')
    # in_segfiles = ' '.join(in_segfiles)

    in_segfiles = os.path.dirname(in_files[0][0])

    cmd = config['python']
    
    cmd_args = [
                merge_segs_script, 
                '--merge_type','outer',
                '--nan_value','NA',
                '--indir', in_segfiles,
                '--key_cols','sample_id',
                '--separator','comma',
                '--type','concatenate',
                '--output', out_file,
                '--regex', '*.segments.csv',
                ]
    
    run_cmd(cmd, cmd_args, max_mem=10, queue=config['queue'], hosts=config['hosts'])


@merge(run_hmmcopy, '{0}/metrics/{1}.corrected_reads.csv'.format(config['out_dir'], config['analysis_id']))
def merge_reads(in_files, out_file):
    make_parent_directory(out_file)
    
    merge_segs_script = os.path.join(bin_dir, 'merge.py')

    # in_readsfiles = glob(os.path.dirname(in_files[0][0]) + '/*.corrected_reads.csv')
    # in_readsfiles = ' '.join(in_readsfiles)

    in_readsfiles = os.path.dirname(in_files[0][0])


    cmd = config['python']
    
    cmd_args = [
                merge_segs_script, 
                '--merge_type','outer',
                '--nan_value','NA',
                '--key_cols','sample_id',
                '--separator','comma',
                '--type','concatenate',
                '--output', out_file,
                '--indir', in_readsfiles,
                '--regex','*.corrected_reads.csv'
                ]
    
    run_cmd(cmd, cmd_args, max_mem=10, queue=config['queue'], hosts=config['hosts'])


@merge(extract_gc_metrics, '{0}/metrics/{1}.gc_metrics.txt'.format(config['out_dir'], config['analysis_id']))
def merge_gc(in_files, out_file):
    make_parent_directory(out_file)
    
    merge_script = os.path.join(bin_dir, 'merge.py')

    # in_gcfiles = glob(os.path.dirname(in_files[0][0]) + '/*.gc_metrics.txt')
    # in_gcfiles = ' '.join(in_gcfiles)

    in_gcfiles = os.path.dirname(in_files[0])


    cmd = config['python']
    
    cmd_args = [
                merge_script, 
                '--merge_type','outer',
                '--nan_value','NA',
                '--indir', in_gcfiles,
                '--key_cols','sample_id',
                '--separator','comma',
                '--type','concatenate',
                '--output', out_file,
                '--regex', '*.gc_metrics.txt',
                ]
    
    run_cmd(cmd, cmd_args, max_mem=10, queue=config['queue'], hosts=config['hosts'])

@follows(extract_metrics_table, extract_quality_metrics)
@files(['{0}/metrics/{1}.quality_metrics.csv'.format(config['out_dir'], config['analysis_id']), r'{0}/metrics/summary/{1}_{2}.metrics.csv'.format(config['out_dir'],library_id, run_id)], '{0}/metrics/summary/allmetrics.csv'.format(config['out_dir']))
# @transform('{0}/metrics/summary/'.format(config['out_dir']), regex('(.*)metrics.csv'), '{0}/metrics/{1}.all_metrics.txt'.format(config['out_dir'], config['analysis_id']))
def merge_metrics(in_files, out_file):
    make_parent_directory(out_file)
    
    cmd = config['python']

    merge_script = os.path.join(bin_dir, 'merge.py')

    cmd_args = [
                merge_script,
                '--merge_type','outer','--nan_value','NA',
                '--key_cols','sample_id','--separator','comma','--type','merge',
                '--output', out_file,
                '--input'
                ]

    cmd_args.extend(in_files)

    
    run_cmd(cmd, cmd_args, queue=config['queue'], hosts=config['hosts'])



@transform(merge_metrics, suffix('.csv'), '.metrics.pdf')
def plot_metrics(metrics_table, out_file):
    make_parent_directory(out_file)
    
    cmd = config['python']
    
    cmd_args = [
                plot_metrics_script,
                metrics_table, 
                out_file,
                '--plot_title','QCpipeline',
                ]
    
    run_cmd(cmd, cmd_args, queue=config['queue'], hosts=config['hosts'],mem=20, max_mem=30)


@follows(merge_metrics, merge_reads)
@files(['{0}/metrics/summary/allmetrics.csv'.format(config['out_dir']), '{0}/metrics/{1}.corrected_reads.csv'.format(config['out_dir'], config['analysis_id'])], '{0}/metrics/summary/heatmap.pdf'.format(config['out_dir']))
def plot_heatmap(inputs, out_file):
    make_parent_directory(out_file)
    
    plot_heatmap_script = os.path.join(bin_dir, 'plot_heatmap.py')

    metrics_table, reads_all = inputs

    cmd = config['python']
    
    cmd_args = [
                plot_heatmap_script,
                '--metrics',metrics_table,
                '--input', reads_all,
                '--separator','comma', 
                '--plot_title','QCpipeline',
                '--column_name','integer_copy_number',
                '--output', out_file
                ]
    
    run_cmd(cmd, cmd_args, queue=config['queue'], hosts=config['hosts'],mem=20, max_mem=30)



@follows(merge_reads, merge_metrics, merge_segs)
@files(['{0}/metrics/{1}.corrected_reads.csv'.format(config['out_dir'], config['analysis_id']),
         '{0}/metrics/{1}.segments.csv'.format(config['out_dir'], config['analysis_id']),
         '{0}/metrics/summary/allmetrics.csv'.format(config['out_dir'])],
        ['{0}/metrics/summary/reads.pdf'.format(config['out_dir']),
        '{0}/metrics/summary/bias.pdf'.format(config['out_dir']),
        '{0}/metrics/summary/segs.pdf'.format(config['out_dir'])])
def plot_hmmcopy(inputs, outputs):
    make_parent_directory(outputs[0])

    plot_hmm_script = os.path.join(bin_dir, 'plot_hmmcopy.py')

    all_reads, segs_all, metrics_table = inputs

    reads_out, bias_output, segs_output = outputs

    cmd = config['python']
    
    cmd_args = [
                plot_hmm_script,
                '--corrected_reads',all_reads,
                '--segments', segs_all,
                '--reference', config['ref_genome'],
                '--quality_metrics', metrics_table,
                '--reads_output', reads_out,
                '--bias_output', bias_output,
                '--segs_output', segs_output, 
                '--plot_title','QCpipeline',
                '--num_states','7'
                ]

    run_cmd(cmd, cmd_args, queue=config['queue'], hosts=config['hosts'],mem=20, max_mem=30)



@follows(plot_metrics, plot_heatmap, plot_hmmcopy)
def end():
    pass

#=======================================================================================================================
# Helper functions
#=======================================================================================================================
def get_wig_chromosomes(file_name):
    chromosomes = []
    
    fh = open(file_name)
    
    for line in fh:
        if 'chrom' not in line:
            continue
        
        line = line.strip()
        
        parts = line.split()
        
        chrom = parts[1]
        
        _, chrom = chrom.split('=')
        
        chromosomes.append(chrom)
    
    fh.close()
    
    return chromosomes

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
