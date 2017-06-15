'''
Extract quality metrics from single cell corrected read counts.
@author: Adi Steif
'''

from __future__ import division

import argparse
import os
import warnings
import numpy as np
import pandas as pd

from statsmodels.robust.scale import stand_mad
from statsmodels.tsa.stattools import acf

#=======================================================================================================================
# Read Command Line Input
#=======================================================================================================================
parser = argparse.ArgumentParser()

parser.add_argument('--hmmcopy_params',
                    help='''Path to HMMcopy result directory with corrected read output .csv files.''')

parser.add_argument('--hmmcopy_corrected_reads',
                    help='''Path to HMMcopy result directory with corrected read output .csv files.''')

parser.add_argument('--hmmcopy_segments',
                    help='''Path to HMMcopy result directory with segments output .csv files.''')

parser.add_argument('--out_file',
                    help='''Path to where output table will be written in .csv format.''')

parser.add_argument('--sample_info',
                    help='''Path to csv file with the sample id information.''')

parser.add_argument('--sample_id',
                    help='''sample id information.''')

args = parser.parse_args()

#=======================================================================================================================
# Functions
#=======================================================================================================================

def compute_total_reads(df):
    total_reads = sum(df['reads'])
    
    return total_reads

def compute_total_reads_hmmcopy(df):
    total_reads_hmmcopy = sum(df['reads'][~df['copy'].isnull()])
    
    return total_reads_hmmcopy

def compute_chr_mad_hmmcopy(df, chr):
    df_hmmcopy = df[~df['copy'].isnull()]
    
    df_chr = df_hmmcopy[df_hmmcopy['chr']==str(chr)]
    
    if len(df_chr) > 0:
        mad_chr = stand_mad(df_chr['copy'], c=1)
    
    else:
        mad_chr = float('NaN')
    
    return mad_chr

def compute_mad_hmmcopy(df):
    df_hmmcopy = df[~df['copy'].isnull()]
    
    if len(df_hmmcopy) > 0:
        mad_hmmcopy = stand_mad(df_hmmcopy['copy'], c=1)
    
    else:
        mad_hmmcopy = float('NaN')
    
    return mad_hmmcopy

def compute_mad_neutral_state(df):
    df_neutral_state = df[(df['state']==3) & (~df['copy'].isnull())]
    
    if len(df_neutral_state) > 0:
        mad_neutral_state = stand_mad(df_neutral_state['copy'], c=1)
    
    else:
        mad_neutral_state = float('NaN')
    
    return mad_neutral_state

def compute_mad_autosomes(df):
    df_autosomes = df[(df['chr']!='X') & (df['chr']!='Y')]
    
    if len(df_autosomes) > 0:
        mad_autosomes = stand_mad(df_autosomes['copy'], c=1)
    
    else:
        mad_autosomes = float('NaN')
    
    return mad_autosomes

def compute_cv_hmmcopy(df):
    df_hmmcopy = df[~df['copy'].isnull()]
    
    if len(df_hmmcopy) > 0:
        df_mean = np.mean(df_hmmcopy['copy'])
        
        df_sd = np.std(df_hmmcopy['copy'])
        
        cv_hmmcopy = df_sd / df_mean
    
    else: 
        cv_hmmcopy = float('NaN')
    
    return cv_hmmcopy

def compute_cv_neutral_state(df):
    df_neutral_state = df[(df['state']==3) & (~df['copy'].isnull())]
    
    if len(df_neutral_state) > 0:
        df_mean = np.mean(df_neutral_state['copy'])
        
        df_sd = np.std(df_neutral_state['copy'])
        
        cv_neutral_state = df_sd / df_mean
    
    else: 
        cv_neutral_state = float('NaN')
    
    return cv_neutral_state

def compute_autocorrelation_hmmcopy(df):
    df_hmmcopy = df[~df['copy'].isnull()]
    
    if len(df_hmmcopy) > 0:
        acf_hmmcopy = acf(df_hmmcopy['copy'], nlags=1)
        
        ac_hmmcopy = acf_hmmcopy[1]
    
    else: 
        ac_hmmcopy = float('NaN')
    
    return ac_hmmcopy

def compute_mean_median_std_hmmcopy_reads_per_bin(df):
    df_hmmcopy = df[~df['copy'].isnull()]
    
    mean_hmmcopy_reads_per_bin = np.mean(df_hmmcopy['reads'])
    
    median_hmmcopy_reads_per_bin = np.median(df_hmmcopy['reads'])
    
    std_hmmcopy_reads_per_bin = np.std(df_hmmcopy['reads'])
    
    return mean_hmmcopy_reads_per_bin, median_hmmcopy_reads_per_bin, std_hmmcopy_reads_per_bin

def compute_num_empty_bins(df):
    df_hmmcopy = df[~df['copy'].isnull()]
    
    empty_bins_hmmcopy = len(df_hmmcopy[df_hmmcopy['reads']==0])
    
    empty_bins_hmmcopy_chrY = len(df_hmmcopy[(df_hmmcopy['reads']==0) & (df_hmmcopy['chr']=='Y')])
    
    return empty_bins_hmmcopy, empty_bins_hmmcopy_chrY

def median_of_bin_residuals_from_segment_median(df, df_seg):
    '''
    MBRSM: This measures "dispersion", but not "integerness".
    '''
    residuals = []
    
    for seg_index in range(len(df_seg)):
        seg_chr = df_seg['chr'][seg_index]
        
        seg_start = df_seg['start'][seg_index]
        
        seg_end = df_seg['end'][seg_index]
        
        seg_median = df_seg['integer_median'][seg_index]
        
        df_seg_bins = df[(df['chr']==seg_chr) & 
                         (df['start'] >= seg_start) & 
                         (df['end'] <= seg_end)]
        
        df_seg_bins_hmmcopy = df_seg_bins[~df_seg_bins['copy'].isnull()]
        
        seg_residuals = [np.abs(x - seg_median) for x in df_seg_bins_hmmcopy['integer_copy_scale']]
        
        residuals.extend(seg_residuals)
    
    median_bin_residuals_median = np.median(residuals)
    
    return median_bin_residuals_median

def median_of_segment_residuals_from_segment_integer(df_seg):
    '''
    MSRSI: This measures "integerness", but not "dispersion".
    '''
    residuals = np.abs(df_seg['integer_median']-df_seg['integer_copy_number'])
    
    median_seg_residuals_integer = np.median(residuals)
    
    return median_seg_residuals_integer

def median_of_bin_residuals_from_segment_integer(df):
    '''
    MBRSI: This measures both "dispersion" and "integerness".
    '''
    df_hmmcopy = df[~df['copy'].isnull()]
    
    residuals = np.abs(df_hmmcopy['integer_copy_scale']-df_hmmcopy['integer_copy_number'])
    
    median_bin_residuals_integer = np.median(residuals)
    
    return median_bin_residuals_integer

def compute_quality_metrics(df, df_seg, sample_id):    
    if 'copy' in df.columns:
        #total_reads = compute_total_reads(df)
        
        total_reads_hmmcopy = compute_total_reads_hmmcopy(df)
        
        mad_chr19 = compute_chr_mad_hmmcopy(df, '19')
        
        mad_hmmcopy = compute_mad_hmmcopy(df)
        
        mad_neutral_state = compute_mad_neutral_state(df)
        
        mad_autosomes = compute_mad_autosomes(df)
        
        cv_hmmcopy = compute_cv_hmmcopy(df)
        
        cv_neutral_state = compute_cv_neutral_state(df)
        
        autocorrelation_hmmcopy = compute_autocorrelation_hmmcopy(df)
        
        mean_hmmcopy_reads_per_bin, median_hmmcopy_reads_per_bin, std_hmmcopy_reads_per_bin = compute_mean_median_std_hmmcopy_reads_per_bin(df)
        
        empty_bins_hmmcopy, empty_bins_hmmcopy_chrY = compute_num_empty_bins(df)
        
        MBRSI_dispersion_non_integerness = median_of_bin_residuals_from_segment_integer(df)
        
        if df_seg is not None:
            MBRSM_dispersion = median_of_bin_residuals_from_segment_median(df, df_seg)
            
            MSRSI_non_integerness = median_of_segment_residuals_from_segment_integer(df_seg)
        
        else:
            MBRSM_dispersion = float('NaN')
            
            MSRSI_non_integerness = float('NaN')
        
        metrics = pd.Series({'sample_id': sample_id, 
#                             'total_reads': total_reads, 
                             'total_reads_hmmcopy': total_reads_hmmcopy,
                             'mad_chr19': mad_chr19, 
                             'mad_hmmcopy': mad_hmmcopy, 
                             'mad_neutral_state': mad_neutral_state, 
                             'mad_autosomes': mad_autosomes, 
                             'cv_hmmcopy': cv_hmmcopy, 
                             'cv_neutral_state': cv_neutral_state, 
                             'autocorrelation_hmmcopy': autocorrelation_hmmcopy, 
                             'mean_hmmcopy_reads_per_bin': mean_hmmcopy_reads_per_bin, 
                             'median_hmmcopy_reads_per_bin': median_hmmcopy_reads_per_bin, 
                             'std_hmmcopy_reads_per_bin': std_hmmcopy_reads_per_bin, 
                             'empty_bins_hmmcopy': empty_bins_hmmcopy, 
                             'empty_bins_hmmcopy_chrY': empty_bins_hmmcopy_chrY, 
                             'MBRSI_dispersion_non_integerness': MBRSI_dispersion_non_integerness,
                             'MBRSM_dispersion': MBRSM_dispersion,
                             'MSRSI_non_integerness': MSRSI_non_integerness})
    
    else: 
        metrics = pd.Series({'sample_id': sample_id, 
 #                            'total_reads': float('NaN'), 
                             'total_reads_hmmcopy': float('NaN'),
                             'mad_chr19': float('NaN'), 
                             'mad_hmmcopy': float('NaN'), 
                             'mad_neutral_state': float('NaN'), 
                             'mad_autosomes': float('NaN'), 
                             'cv_hmmcopy': float('NaN'), 
                             'cv_neutral_state': float('NaN'), 
                             'autocorrelation_hmmcopy': float('NaN'), 
                             'mean_hmmcopy_reads_per_bin': float('NaN'), 
                             'median_hmmcopy_reads_per_bin': float('NaN'), 
                             'std_hmmcopy_reads_per_bin': float('NaN'), 
                             'empty_bins_hmmcopy': float('NaN'), 
                             'empty_bins_hmmcopy_chrY': float('NaN'), 
                             'MBRSI_dispersion_non_integerness': float('NaN'),
                             'MBRSM_dispersion': float('NaN'),
                             'MSRSI_non_integerness': float('NaN')})
    
    return metrics


def extract_sample_info(sample_info, sample_id, infile):
    """
    get info
    """

    if not sample_info and sample_id:
        return sample_id

    if not sample_info and not sample_id:
        return os.path.basename(args.hmmcopy_corrected_reads).split('.')[0]

    sample_info = open(sample_info)
    header = sample_info.readline()
    sampdata = sample_info.readline()

    assert sample_info.readline() == ''

    sampdata = sampdata.strip().split(',')
    samp = sampdata[0]

    return samp

#=======================================================================================================================
# Run script
#=======================================================================================================================
'''
args.hmmcopy_results_dir = '/share/scratch/asteif_temp/single_cell_indexing/hmmcopy/PX0286_merge_C7H39ANXX_7_C7PJ9ANXX_7/SA501X3F-1/2016-02-03_bin_200000_e_995_s_35_kappa_diploid/results'
args.out_file = '/share/lustre/asteif/projects/single_cell_indexing/test/hmmcopy_heatmap/PX0286_merge_C7H39ANXX_7_C7PJ9ANXX_7.quality_metrics.csv'

Note 1: The metric 'mad_hmmcopy', which is the MAD for all non-missing positions, 
does not work well for low-coverage samples. Those samples have so little coverage
that the 'mad_hmmcopy' is exactly 0. Consider removing this metric.

Note 2: PCA is sensitive to the relative scaling of the original variables, 
so if number of reads is included, it will take up most of the variance
should apply mean centering and normalizing for each attribute

- Total number of reads
- Total reads in HMMcopy bins (not NA in 'copy')
- MAD (median absolute deviation) for HMMcopy bins
- Neutral MAD (HMMcopy non-NA bins in state 3)
- CV (standard deviation divided by the mean) for HMMcopy bins
- Neutral CV (HMMcopy non-NA bins in state 3)
- Autocorrelation coefficient (lag of one bin)
- Log likelihood
- Mean number of reads per bin for HMMcopy bins
- Median number of reads per bin for HMMcopy bins
- Standard deviation in number of reads per bin for HMMcopy bins
- Number of empty bins for HMMcopy bins
- Number of empty bins for HMMcopy bins in chrY
'''

def main():
    df_metrics = pd.DataFrame()
    
    df_loglik = pd.Series()

    sample_id = extract_sample_info(args.sample_info, args.sample_id, args.hmmcopy_corrected_reads)

    corrected_data = pd.read_csv(args.hmmcopy_corrected_reads)

    if os.stat(args.hmmcopy_segments).st_size != 0:
        segments_data = pd.read_csv(args.hmmcopy_segments)
    else:
        segments_data = None


    metrics = compute_quality_metrics(corrected_data, segments_data, sample_id)

    df_metrics = pd.concat([df_metrics, pd.DataFrame(metrics).transpose()])


    if os.stat(args.hmmcopy_params).st_size != 0:
        param_data = pd.read_csv(args.hmmcopy_params)

        loglik = float(param_data.ix[param_data['parameter']=='loglik', 'final'])

        df_loglik = df_loglik.append(pd.Series({sample_id: loglik}))

    else:
        df_loglik = df_loglik.append(pd.Series({sample_id: np.nan}))

    df_loglik = pd.DataFrame(df_loglik.reset_index())
    
    df_loglik.columns = ['sample_id', 'log_likelihood']
    
    df_metrics = df_metrics.merge(df_loglik)
    
    df_metrics.sort_values('sample_id', inplace=True)
    
    df_metrics.reset_index(inplace=True, drop=True)
    
    df_metrics.to_csv(args.out_file, index=False)

if __name__ == '__main__':
    main()

