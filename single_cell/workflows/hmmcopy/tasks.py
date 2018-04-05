'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner
from scripts import ExtractHmmMetrics
from scripts import GenerateCNMatrix
from scripts import FilterHmmData
from scripts import ConvertCSVToSEG
from scripts import ReadCounter
from scripts import CorrectReadCount

from single_cell.utils import csvutils
from single_cell.utils import pdfutils
from single_cell.utils import helpers

scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
run_hmmcopy_rscript = os.path.join(scripts_directory, 'hmmcopy.R')


def run_hmmcopy(
    bam_file,
    bai_file,
    corrected_reads_filename,
    segments_filename,
    parameters_filename,
    posterior_marginals_filename,
    hmm_metrics,
    cell_id,
    config,
    hmmparams,
    tempdir,):


    run_readcount_rscript = os.path.join(scripts_directory, 'correct_read_count.R')

    #generate wig file for hmmcopy
    os.makedirs(tempdir)
    readcount_wig = os.path.join(tempdir, 'readcounter.wig')

    rc = ReadCounter(bam_file, readcount_wig, hmmparams['bin_size'], config['chromosomes'],
                     hmmparams['min_mqual'], excluded=hmmparams['exclude_list'])
    rc.main()

    correct_reads_out = os.path.join(tempdir, 'corrected_reads.csv')
    if hmmparams["smoothing_function"] == 'loess':
        cmd=['Rscript', run_readcount_rscript,
             readcount_wig,
             hmmparams['gc_wig_file'],
             hmmparams['map_wig_file'],
             correct_reads_out
             ]
        pypeliner.commandline.execute(*cmd)
    elif hmmparams["smoothing_function"] == 'modal':
        CorrectReadCount(hmmparams["gc_wig_file"],
                         hmmparams['map_wig_file'],
                         readcount_wig,
                         correct_reads_out,
                         mappability=hmmparams['map_cutoff']).main()
    else:
        raise Exception("smoothing function %s not supported. pipeline supports loess and modal" %hmmparams["smoothing_function"])

    #run hmmcopy
    cmd = ['Rscript', run_hmmcopy_rscript,
        '--corrected_data=' + correct_reads_out,
        '--reads_output=' + corrected_reads_filename,
        '--segs_output=' + segments_filename,
        '--params_output=' + parameters_filename,
        '--post_marginals_output=' + posterior_marginals_filename,
        '--sample_id=' + cell_id]

    if hmmparams['map_cutoff']:
        cmd.append('--map_cutoff=' + str(hmmparams['map_cutoff']))

    if hmmparams['num_states']:
        cmd.append('--num_states=' + str(hmmparams['num_states']))

    if hmmparams['mu']:
        cmd.append('--param_mu=' + str(hmmparams['mu']))

    if hmmparams['m']:
        cmd.append('--param_m=' + str(hmmparams['m']))

    if hmmparams['kappa']:
        cmd.append('--param_k=' + str(hmmparams['kappa']))

    if hmmparams['e']:
        cmd.append('--param_e=' + str(hmmparams['e']))

    if hmmparams['g']:
        cmd.append('--param_g=' + str(hmmparams['g']))

    if hmmparams['s']:
        cmd.append('--param_s=' + str(hmmparams['s']))

    if hmmparams['strength']:
        cmd.append('--param_str=' + str(hmmparams['strength']))

    if hmmparams['nu']:
        cmd.append('--param_nu=' + str(hmmparams['nu']))

    if hmmparams['eta']:
        cmd.append('--param_eta=' + str(hmmparams['eta']))

    if hmmparams['lambda']:
        cmd.append('--param_l=' + str(hmmparams['lambda']))

    if hmmparams["auto_ploidy"]:
        cmd.append('--auto_ploidy')

    pypeliner.commandline.execute(*cmd)

    # generate the metrics file for hmmcopy
    metrics = ExtractHmmMetrics(parameters_filename, corrected_reads_filename,
                                segments_filename, hmm_metrics, cell_id)
    metrics.main()


def merge_files(reads, segs, hmm_metrics, hmm_params, hmm_posteriors,
                merged_segs, merged_reads, merged_hmm_metrics,
                merged_hmm_params, merged_hmm_posteriors, cn_matrix, temp,
                mad_thres,merged_reads_filt, merged_segs_filt, igv_segs):

    csvutils.concatenate_csv_lowmem(reads, merged_reads)
    csvutils.concatenate_csv_lowmem(segs, merged_segs)
    csvutils.concatenate_csv_lowmem(hmm_metrics, merged_hmm_metrics)

    csvutils.concatenate_csv(hmm_params, merged_hmm_params)
    csvutils.concatenate_csv_lowmem(hmm_posteriors, merged_hmm_posteriors)

    generate_cn_matrix(reads, cn_matrix, temp)
    
    filter_hmm_data(merged_hmm_metrics, merged_segs, merged_reads, mad_thres,
                    merged_reads_filt, merged_segs_filt)

    convert_csv_to_seg(merged_segs_filt, merged_reads_filt, igv_segs)


def generate_cn_matrix(infiles, output, tempdir):
    """
    generate cn matrix per cell
    merge all matrices and dump to output
    """
    
    if not os.path.exists(tempdir):
        os.mkdir(tempdir)
    
    matrix_files = []
    for cell_id, infile in infiles.iteritems():
        outfile= os.path.join(tempdir, cell_id+'_matrix.csv')
        gen_gc = GenerateCNMatrix(infile, outfile, ',', 'integer_copy_number',
                                  cell_id, 'hmmcopy_corrected_reads')
        gen_gc.main()
        matrix_files.append(outfile)

    csvutils.merge_csv(matrix_files, output, 'outer', 'chr,start,end,width')
    
    
def filter_hmm_data(quality_metrics, segments, reads, mad_threshold,
                    reads_out, segments_out):
    filter_hmm = FilterHmmData(quality_metrics, segments, reads,
                               mad_threshold, reads_out, segments_out)
    filter_hmm.main()


def convert_csv_to_seg(filtered_segs, filtered_reads, output_seg):
    converter = ConvertCSVToSEG(filtered_segs, filtered_reads, output_seg)
    converter.main()
