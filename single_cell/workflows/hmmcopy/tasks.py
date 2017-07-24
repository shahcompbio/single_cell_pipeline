'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner

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
        config['Rscript'], run_hmmcopy_rscript,
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
