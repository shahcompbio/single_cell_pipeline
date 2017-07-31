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


def run_hmmcopy(bam_file,
    corrected_reads_filename,
    segments_filename,
    parameters_filename,
    posterior_marginals_filename,
    hmm_metrics,
    sample_id,
    config,
    tempdir,):

    #generate wig file for hmmcopy
    os.makedirs(tempdir)
    readcount_wig = os.path.join(tempdir, 'readcounter.wig')
    pypeliner.commandline.execute(
            config['readcounter'],
            '-w', str(config['bin_size']),
            '-q', str(config['min_mqual']),
            '-c', ','.join(config['chromosomes']),
            bam_file,
            '>',
            readcount_wig,
        )


    #run hmmcopy
    cmd = [config['Rscript'], run_hmmcopy_rscript,
        '--tumour_file=' + readcount_wig,
        '--gc_file=' + config['gc_wig_file'],
        '--map_file=' + config['map_wig_file'],
        '--reads_output=' + corrected_reads_filename,
        '--segs_output=' + segments_filename,
        '--params_output=' + parameters_filename,
        '--post_marginals_output=' + posterior_marginals_filename,
        '--sample_id=' + sample_id]

    if config['parameters']['map_cutoff']:
        cmd.append('--map_cutoff=' + str(config['parameters']['map_cutoff']))

    if config['parameters']['num_states']:
        cmd.append('--num_states=' + str(config['parameters']['num_states']))

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

    # generate the metrics file for hmmcopy
    metrics = ExtractHmmMetrics(parameters_filename, corrected_reads_filename,
                                segments_filename, hmm_metrics, sample_id)
    metrics.main()


def generate_cn_matrix(infile, output, sep, colname, sample_id, typ):
    gen_gc = GenerateCNMatrix(infile, output, sep, colname, sample_id, typ)
    gen_gc.main()
