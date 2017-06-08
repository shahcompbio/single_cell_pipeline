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

parser.add_argument('--hmmcopy_results_dir',
                    help='''Path to HMMcopy result directory with corrected read output .csv files.''')

parser.add_argument('--out_file',
                    help='''Path to where output table will be written in .csv format.''')

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

def compute_quality_metrics(df, sample_id):
    if 'copy' in df.columns:
        total_reads = compute_total_reads(df)
        
        total_reads_hmmcopy = compute_total_reads_hmmcopy(df)
        
        mad_chr19 = compute_chr_mad_hmmcopy(df, '19')
        
        mad_hmmcopy = compute_mad_hmmcopy(df)
        
        mad_neutral_state = compute_mad_neutral_state(df)
        
        cv_hmmcopy = compute_cv_hmmcopy(df)
        
        cv_neutral_state = compute_cv_neutral_state(df)
        
        autocorrelation_hmmcopy = compute_autocorrelation_hmmcopy(df)
        
        metrics = pd.Series({'sample_id': sample_id, 
                             'total_reads': total_reads, 
                             'total_reads_hmmcopy': total_reads_hmmcopy,
                             'mad_chr19': mad_chr19, 
                             'mad_hmmcopy': mad_hmmcopy, 
                             'mad_neutral_state': mad_neutral_state, 
                             'cv_hmmcopy': cv_hmmcopy, 
                             'cv_neutral_state': cv_neutral_state, 
                             'autocorrelation_hmmcopy': autocorrelation_hmmcopy})
    
    else: 
        metrics = pd.Series({'sample_id': sample_id, 
                             'total_reads': float('NaN'), 
                             'total_reads_hmmcopy': float('NaN'),
                             'mad_chr19': float('NaN'), 
                             'mad_hmmcopy': float('NaN'), 
                             'mad_neutral_state': float('NaN'), 
                             'cv_hmmcopy': float('NaN'), 
                             'cv_neutral_state': float('NaN'), 
                             'autocorrelation_hmmcopy': float('NaN')})
    
    return metrics

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

'''

def main():
    df_metrics = pd.DataFrame()
    
    df_loglik = pd.Series()
    
    for file in os.listdir(args.hmmcopy_results_dir):
        if file.endswith('.corrected_reads.csv'):
            sample_id = os.path.basename(file).split('.')[0]
            
            corrected_data = pd.read_csv(os.path.join(args.hmmcopy_results_dir, file))
            
            metrics = compute_quality_metrics(corrected_data, sample_id)
            
            df_metrics = pd.concat([df_metrics, pd.DataFrame(metrics).transpose()])
            
        elif file.endswith('.parameters.csv'):
            sample_id = os.path.basename(file).split('.')[0]
            
            if os.stat(os.path.join(args.hmmcopy_results_dir, file)).st_size != 0:
                param_data = pd.read_csv(os.path.join(args.hmmcopy_results_dir, file))
                
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
