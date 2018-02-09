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
from single_cell.utils import concatenate_csv

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
    pypeliner.commandline.execute(
            'readCounter',
            '-w', str(hmmparams['bin_size']),
            '-q', str(hmmparams['min_mqual']),
            '-c', ','.join(config['chromosomes']),
            bam_file,
            '>',
            readcount_wig,
        )

    #run hmmcopy
    cmd = ['Rscript', run_hmmcopy_rscript,
        '--tumour_file=' + readcount_wig,
        '--gc_file=' + hmmparams['gc_wig_file'],
        '--map_file=' + hmmparams['map_wig_file'],
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


def concatenate_csv_pandas(in_filenames, out_filename, nan_val='NA'):
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



def merge_files(reads, segs, hmm_metrics, hmm_params, hmm_posteriors,
                merged_segs, merged_reads, merged_hmm_metrics,
                merged_hmm_params, merged_hmm_posteriors, cn_matrix, temp,
                mad_thres,merged_reads_filt, merged_segs_filt, igv_segs):

    concatenate_csv(reads, merged_reads)
    concatenate_csv(segs, merged_segs)
    concatenate_csv(hmm_metrics, merged_hmm_metrics)

    concatenate_csv_pandas(hmm_params, merged_hmm_params)
    concatenate_csv(hmm_posteriors, merged_hmm_posteriors)

    generate_cn_matrix(reads, cn_matrix, temp)
    
    filter_hmm_data(merged_hmm_metrics, merged_segs, merged_reads, mad_thres,
                    merged_reads_filt, merged_segs_filt)

    convert_csv_to_seg(merged_segs_filt, merged_reads_filt, igv_segs)


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

def plot_hmmcopy(reads, segments, metrics, sample_info, ref_genome, reads_out, segs_out,
                 bias_out, sample_id, num_states=7, plot_title=None,
                 mad_threshold=None, annotation_cols=None):
    plot = GenHmmPlots(reads, segments, metrics, sample_info, ref_genome, reads_out, segs_out,
                       bias_out, sample_id, num_states=num_states, plot_title=plot_title,
                       mad_threshold=mad_threshold, annotation_cols=annotation_cols)
    plot.main()

def merge_pdf(in_filenames, out_filename, metrics, mad_threshold):

    metrics = pd.read_csv(metrics, sep=',')

    data = {}
    for cell, cond in zip(metrics['cell_id'], metrics['experimental_condition']):
        if cond in data:
            data[cond].append(cell)
        else:
            data[cond] = [cell]
    samples = [x for _,v in data.iteritems() for x in sorted(v)]
    
    for in_files, out_file in zip(in_filenames, out_filename):

        outdir = os.path.dirname(out_file)
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        merger = PdfFileMerger()

        for samp in samples:
            infile = in_files[samp]
            #filter by mad if mad_threshold is specified
            if mad_threshold:
                mad = metrics[metrics['cell_id'] == samp]['mad_neutral_state'].iloc[0]
                if mad>mad_threshold:
                    continue
            
            merger.append(open(infile, 'rb'))
        if not os.path.exists(os.path.dirname(out_file)):
            os.makedirs(os.path.dirname(out_file))
        with open(out_file, 'wb') as fout:
            merger.write(fout)

def convert_csv_to_seg(filtered_segs, filtered_reads, output_seg):
    converter = ConvertCSVToSEG(filtered_segs, filtered_reads, output_seg)
    converter.main()
