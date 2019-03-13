'''
Extract quality metrics from single cell corrected read counts.
@author: Adi Steif
'''

from __future__ import division

import os
import numpy as np
import pandas as pd
import logging

from statsmodels.robust.scale import stand_mad
from statsmodels.tsa.stattools import acf

#=========================================================================
# Functions
#=========================================================================


class ExtractHmmMetrics(object):

    def __init__(self, params, reads, segments, output, sample_id, table_name = None):
        self.params = params
        self.segments = segments
        self.reads = reads
        self.output = output
        self.sample_id = sample_id
        self.output_format = self.get_output_format(self.output)
        self.table_name = table_name

        if not self.table_name:
            self.table_name = 'hmmcopy/metrics/{}'.format(self.sample_id)

    def get_output_format(self, output):
        _, ext = os.path.splitext(output)

        if ext == ".csv":
            return "csv"
        elif ext == ".gz":
            return "gzip"
        elif ext == ".h5" or ext == ".hdf5":
            return "h5"
        else:
            logging.getLogger("single_cell.hmmcopy.extract_metrics").warn(
                "Couldn't detect output format. extension {}".format(ext))
            return "csv"

    def write_df_to_file(self, df):

        if self.output_format == "csv":
            df.to_csv(self.output, na_rep='NA', index=False)

        elif self.output_format == "gzip":
            df.to_csv(self.output, na_rep='NA', index=False, compression="gzip")
        else:
            df = df.infer_objects()

            with pd.HDFStore(self.output, 'w', complevel=9, complib='blosc') as out_store:
                out_store.put(self.table_name, df, format='table')

    def compute_total_reads(self, df):
        total_reads = sum(df['reads'])

        return total_reads

    def compute_total_reads_hmmcopy(self, df):
        total_reads_hmmcopy = sum(df['reads'][~df['copy'].isnull()])

        return total_reads_hmmcopy

    def compute_chr_mad_hmmcopy(self, df, chrom):
        df_hmmcopy = df[~df['copy'].isnull()]

        df_chr = df_hmmcopy[df_hmmcopy['chr'] == str(chrom)]

        if len(df_chr) > 0:
            mad_chr = stand_mad(df_chr['copy'], c=1)

        else:
            mad_chr = float('NaN')

        return mad_chr

    def compute_mad_hmmcopy(self, df):
        df_hmmcopy = df[~df['copy'].isnull()]

        if len(df_hmmcopy) > 0:
            mad_hmmcopy = stand_mad(df_hmmcopy['copy'], c=1)

        else:
            mad_hmmcopy = float('NaN')

        return mad_hmmcopy

    def compute_mad_neutral_state(self, df):
        df_neutral_state = df[(df['state'] == 3) & (~df['copy'].isnull())]

        if len(df_neutral_state) > 0:
            mad_neutral_state = stand_mad(df_neutral_state['copy'], c=1)

        else:
            mad_neutral_state = float('NaN')

        return mad_neutral_state

    def compute_mad_autosomes(self, df):
        df_autosomes = df[(df['chr'] != 'X') & (df['chr'] != 'Y')]

        if len(df_autosomes) > 0:
            mad_autosomes = stand_mad(df_autosomes['copy'], c=1)

        else:
            mad_autosomes = float('NaN')

        return mad_autosomes

    def compute_cv_hmmcopy(self, df):
        df_hmmcopy = df[~df['copy'].isnull()]

        if len(df_hmmcopy) > 0:
            df_mean = np.mean(df_hmmcopy['copy'])

            df_sd = np.std(df_hmmcopy['copy'])

            cv_hmmcopy = df_sd / df_mean

        else:
            cv_hmmcopy = float('NaN')

        return cv_hmmcopy

    def compute_cv_neutral_state(self, df):
        df_neutral_state = df[(df['state'] == 3) & (~df['copy'].isnull())]

        if len(df_neutral_state) > 0:
            df_mean = np.mean(df_neutral_state['copy'])

            df_sd = np.std(df_neutral_state['copy'])

            cv_neutral_state = df_sd / df_mean

        else:
            cv_neutral_state = float('NaN')

        return cv_neutral_state

    def compute_autocorrelation_hmmcopy(self, df):
        df_hmmcopy = df[~df['copy'].isnull()]

        if len(df_hmmcopy) > 0:
            acf_hmmcopy = acf(df_hmmcopy['copy'], nlags=1)

            ac_hmmcopy = acf_hmmcopy[1]

        else:
            ac_hmmcopy = float('NaN')

        return ac_hmmcopy

    def compute_mean_median_std_hmmcopy_reads_per_bin(self, df):
        df_hmmcopy = df[~df['copy'].isnull()]

        mean_hmmcopy_reads_per_bin = np.mean(df_hmmcopy['reads'])

        median_hmmcopy_reads_per_bin = np.median(df_hmmcopy['reads'])

        std_hmmcopy_reads_per_bin = np.std(df_hmmcopy['reads'])

        return mean_hmmcopy_reads_per_bin, median_hmmcopy_reads_per_bin, std_hmmcopy_reads_per_bin

    def compute_num_empty_bins(self, df):
        df_hmmcopy = df[~df['copy'].isnull()]

        empty_bins_hmmcopy = len(df_hmmcopy[df_hmmcopy['reads'] == 0])

        empty_bins_hmmcopy_chrY = len(
            df_hmmcopy[
                (df_hmmcopy['reads'] == 0) & (
                    df_hmmcopy['chr'] == 'Y')])

        return empty_bins_hmmcopy, empty_bins_hmmcopy_chrY

    def median_of_bin_residuals_from_segment_median(self, df, df_seg):
        '''
        MBRSM: This measures "dispersion", but not "integerness".
        '''
        residuals = []

        for seg_index in range(len(df_seg)):
            seg_chr = df_seg['chr'][seg_index]

            seg_start = df_seg['start'][seg_index]

            seg_end = df_seg['end'][seg_index]

            seg_median = df_seg['integer_median'][seg_index]

            df_seg_bins = df[(df['chr'] == seg_chr) &
                             (df['start'] >= seg_start) &
                             (df['end'] <= seg_end)]

            df_seg_bins_hmmcopy = df_seg_bins[~df_seg_bins['copy'].isnull()]

            seg_residuals = [
                np.abs(
                    x -
                    seg_median) for x in df_seg_bins_hmmcopy['integer_copy_scale']]

            residuals.extend(seg_residuals)

        median_bin_residuals_median = np.median(residuals)

        return float(median_bin_residuals_median)

    def median_of_segment_residuals_from_segment_integer(self, df_seg):
        '''
        MSRSI: This measures "integerness", but not "dispersion".
        '''
        residuals = np.abs(
            df_seg['integer_median'] -
            df_seg['integer_copy_number'])

        if len(residuals):
            median_seg_residuals_integer = np.median(residuals)
        else:
            median_seg_residuals_integer = float('nan')

        return float(median_seg_residuals_integer)

    def median_of_bin_residuals_from_segment_integer(self, df):
        '''
        MBRSI: This measures both "dispersion" and "integerness".
        '''
        df_hmmcopy = df[~df['copy'].isnull()]

        residuals = np.abs(
            df_hmmcopy['integer_copy_scale'] -
            df_hmmcopy['integer_copy_number'])

        median_bin_residuals_integer = np.median(residuals)

        return float(median_bin_residuals_integer)

    def compute_quality_metrics(self, df, df_seg, sample_id):
        if 'copy' in df.columns:
            total_reads = self.compute_total_reads(df)

            total_reads_hmmcopy = self.compute_total_reads_hmmcopy(df)

            mad_chr19 = self.compute_chr_mad_hmmcopy(df, '19')

            mad_hmmcopy = self.compute_mad_hmmcopy(df)

            mad_neutral_state = self.compute_mad_neutral_state(df)

            mad_autosomes = self.compute_mad_autosomes(df)

            cv_hmmcopy = self.compute_cv_hmmcopy(df)

            cv_neutral_state = self.compute_cv_neutral_state(df)

            autocorrelation_hmmcopy = self.compute_autocorrelation_hmmcopy(df)

            mean_hmmcopy_reads_per_bin, median_hmmcopy_reads_per_bin, std_hmmcopy_reads_per_bin = self.compute_mean_median_std_hmmcopy_reads_per_bin(
                df)

            empty_bins_hmmcopy, empty_bins_hmmcopy_chrY = self.compute_num_empty_bins(
                df)

            MBRSI_dispersion_non_integerness = self.median_of_bin_residuals_from_segment_integer(
                df)

            if df_seg is not None:
                MBRSM_dispersion = self.median_of_bin_residuals_from_segment_median(
                    df,
                    df_seg)

                MSRSI_non_integerness = self.median_of_segment_residuals_from_segment_integer(
                    df_seg)

            else:
                MBRSM_dispersion = float('NaN')

                MSRSI_non_integerness = float('NaN')

            metrics = pd.Series({'cell_id': sample_id,
                                 'total_mapped_reads_hmmcopy': total_reads,
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
            metrics = pd.Series({'cell_id': sample_id,
                                 'total_mapped_reads_hmmcopy': float('NaN'),
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

    def main(self):
        df_metrics = pd.DataFrame()

        df_loglik = pd.Series()

        corrected_data = pd.read_csv(self.reads)

        if os.stat(self.segments).st_size != 0:
            segments_data = pd.read_csv(self.segments)
        else:
            segments_data = None

        metrics = self.compute_quality_metrics(
            corrected_data,
            segments_data,
            self.sample_id)

        df_metrics = pd.concat([df_metrics, pd.DataFrame(metrics).transpose()])

        param_data = pd.read_csv(self.params)

        if not param_data.empty:

            loglik = float(
                param_data.ix[
                    param_data['parameter'] == 'loglik',
                    'final'])

            df_loglik = df_loglik.append(pd.Series({self.sample_id: loglik}))
        else:
            df_loglik = df_loglik.append(pd.Series({self.sample_id: np.nan}))

        df_loglik = pd.DataFrame(df_loglik.reset_index())

        df_loglik.columns = ['cell_id', 'log_likelihood']

        df_metrics = df_metrics.merge(df_loglik)

        df_metrics.sort_values('cell_id', inplace=True)

        df_metrics.reset_index(inplace=True, drop=True)

        self.write_df_to_file(df_metrics)
