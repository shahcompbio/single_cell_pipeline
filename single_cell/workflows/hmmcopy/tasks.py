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
    sample_id,
    config,
    hmmparams,
    tempdir,):

    #generate wig file for hmmcopy
    os.makedirs(tempdir)
    readcount_wig = os.path.join(tempdir, 'readcounter.wig')

    rc = ReadCounter(bam_file, readcount_wig, hmmparams['bin_size'], config['chromosomes'],
                     hmmparams['min_mqual'], excluded=hmmparams['exclude_list'])
    rc.main()

    correct_reads_out = os.path.join(tempdir, 'corrected_reads.csv')
    CorrectReadCount(hmmparams["gc_wig_file"],
                     hmmparams['map_wig_file'],
                     readcount_wig,
                     correct_reads_out,
                     mappability=hmmparams['map_cutoff']).main()
    
    #run hmmcopy
    cmd = ['Rscript', run_hmmcopy_rscript,
        '--corrected_data=' + correct_reads_out,
        '--reads_output=' + corrected_reads_filename,
        '--segs_output=' + segments_filename,
        '--params_output=' + parameters_filename,
        '--post_marginals_output=' + posterior_marginals_filename,
        '--sample_id=' + sample_id]

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

    if hmmparams['s']:
        cmd.append('--param_s=' + str(hmmparams['s']))

    pypeliner.commandline.execute(*cmd)

    # generate the metrics file for hmmcopy
    metrics = ExtractHmmMetrics(parameters_filename, corrected_reads_filename,
                                segments_filename, hmm_metrics, sample_id)
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
    for sample_id, infile in infiles.iteritems():
        outfile= os.path.join(tempdir, sample_id+'_matrix.csv')
        gen_gc = GenerateCNMatrix(infile, outfile, ',', 'integer_copy_number',
                                  sample_id, 'hmmcopy_corrected_reads')
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
