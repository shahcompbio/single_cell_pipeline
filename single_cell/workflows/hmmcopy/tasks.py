'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner
from scripts import ExtractHmmMetrics
from scripts import GenerateCNMatrix
from scripts import FilterHmmData
from scripts import GenHmmPlots
from scripts import ConvertCSVToSEG
import pandas as pd
from PyPDF2 import PdfFileMerger


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
            'readCounter',
            '-w', str(config['hmmcopy_params']['bin_size']),
            '-q', str(config['hmmcopy_params']['min_mqual']),
            '-c', ','.join(config['chromosomes']),
            bam_file,
            '>',
            readcount_wig,
        )

    #run hmmcopy
    cmd = ['Rscript', run_hmmcopy_rscript,
        '--tumour_file=' + readcount_wig,
        '--gc_file=' + config['gc_wig_file'],
        '--map_file=' + config['map_wig_file'],
        '--reads_output=' + corrected_reads_filename,
        '--segs_output=' + segments_filename,
        '--params_output=' + parameters_filename,
        '--post_marginals_output=' + posterior_marginals_filename,
        '--sample_id=' + sample_id]

    if config['hmmcopy_params']['map_cutoff']:
        cmd.append('--map_cutoff=' + str(config['hmmcopy_params']['map_cutoff']))

    if config['hmmcopy_params']['num_states']:
        cmd.append('--num_states=' + str(config['hmmcopy_params']['num_states']))

    if config['hmmcopy_params']['mu']:
        cmd.append('--param_mu=' + str(config['hmmcopy_params']['mu']))

    if config['hmmcopy_params']['m']:
        cmd.append('--param_m=' + str(config['hmmcopy_params']['m']))

    if config['hmmcopy_params']['kappa']:
        cmd.append('--param_k=' + str(config['hmmcopy_params']['kappa']))

    if config['hmmcopy_params']['e']:
        cmd.append('--param_e=' + str(config['hmmcopy_params']['e']))

    if config['hmmcopy_params']['s']:
        cmd.append('--param_s=' + str(config['hmmcopy_params']['s']))

    pypeliner.commandline.execute(*cmd)

    # generate the metrics file for hmmcopy
    metrics = ExtractHmmMetrics(parameters_filename, corrected_reads_filename,
                                segments_filename, hmm_metrics, sample_id)
    metrics.main()


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

def merge_csv(in_filenames, out_filename, how, on, nan_val='NA'):
    data = []

    if isinstance(in_filenames, dict):
        in_filenames = in_filenames.values()
    
    for in_filename in in_filenames:
        with open(in_filename) as f:
            first_line = f.readline()
            if len(first_line) == 0:
                continue
        data.append(pd.read_csv(in_filename, dtype=str))

    data = merge_frames(data, how, on)
    data = data.fillna(nan_val)
    data.to_csv(out_filename, index=False)

def concatenate_csv(in_filenames, out_filename, nan_val='NA'):
    data = []
    for _, in_filename in in_filenames.iteritems():
        with open(in_filename) as f:
            first_line = f.readline()
            if len(first_line) == 0:
                continue
        data.append(pd.read_csv(in_filename, dtype=str))
    data = pd.concat(data, ignore_index=True)
    data = data.fillna(nan_val)
    data.to_csv(out_filename, index=False)



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

    merge_csv(matrix_files, output, 'outer', 'chr,start,end,width')
    
    
def filter_hmm_data(quality_metrics, segments, reads, mad_threshold,
                    reads_out, segments_out):
    filter_hmm = FilterHmmData(quality_metrics, segments, reads,
                               mad_threshold, reads_out, segments_out)
    filter_hmm.main()

def plot_hmmcopy(reads, segments, metrics, ref_genome, reads_out, segs_out,
                 bias_out, sample_id, num_states=7, plot_title=None, mad_threshold=None):
    plot = GenHmmPlots(reads, segments, metrics, ref_genome, reads_out, segs_out,
                       bias_out, sample_id, num_states=num_states, plot_title=plot_title,
                       mad_threshold=mad_threshold)
    plot.main()


def merge_pdf(in_filenames, out_filename, metrics, mad_threshold):

    metrics = pd.read_csv(metrics, sep=',')
    
    for in_files, out_file in zip(in_filenames, out_filename):

        merger = PdfFileMerger()



        for samp, infile in in_files.iteritems():
            #filter by mad if mad_threshold is specified
            if mad_threshold:
                mad = metrics[metrics['cell_id'] == samp]['mad_neutral_state'].iloc[0]
                if mad>mad_threshold:
                    continue
            
            merger.append(open(infile, 'rb'))
        
        with open(out_file, 'wb') as fout:
            merger.write(fout)

def convert_csv_to_seg(filtered_segs, filtered_reads, output_seg):
    converter = ConvertCSVToSEG(filtered_segs, filtered_reads, output_seg)
    converter.main()
