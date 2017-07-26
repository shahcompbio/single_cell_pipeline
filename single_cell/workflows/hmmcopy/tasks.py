'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner
from scripts import ExtractHmmMetrics
from scripts import GenerateCNMatrix


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


    cmd = [config['Rscript'], run_hmmcopy_rscript,
        '--tumour_file=' + readcount_wig_filename,
        '--gc_file=' + config['gc_wig_file'],
        '--map_file=' + config['map_wig_file'],
        '--reads_output=' + corrected_reads_filename,
        '--segs_output=' + segments_filename,
        '--params_output=' + parameters_filename,
        '--post_marginals_output=' + posterior_marginals_filename,
        '--sample_id=' + sample_id]

    if config['map_cutoff']:
        cmd.append('--map_cutoff=' + str(config['map_cutoff']))

    if config['num_states']:
        cmd.append('--num_states=' + str(config['num_states']))

    if config['parameters']['mu']:
        cmd.append('--param_mu=' + str(config['parameters']['mu']))

    if config['parameters']['m']:
        cmd.append('--param_m=' + str(config['parameters']['m']))

    if config['parameters']['kappa']:
        cmd.append('--param_k=' + str(config['parameters']['kappa']))

    if config['parameters']['e']:
        cmd.append('--param_e=' + str(config['parameters']['e']))

    if config['parameters']['s']:
        cmd.append('--param_s=' + str(config['parameters']['s']))

    pypeliner.commandline.execute(*cmd)

def extract_hmm_metrics(params, reads, segments, output, sample_id):
    metrics = ExtractHmmMetrics(params, reads, segments, output, sample_id)
    metrics.main()


def collect_cn_metrics(infile, output, sep, colname, sample_id, typ):
    gen_gc = GenerateCNMatrix(infile, output, sep, colname, sample_id, typ)
    gen_gc.main()