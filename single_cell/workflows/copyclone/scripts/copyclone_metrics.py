'''
Created on Mar 28, 2018

@author: dgrewal
'''
import numpy as np
from statsmodels.robust.scale import stand_mad
from statsmodels.tsa.stattools import acf
from math import log


def compute_halfiness(df, df_seg):
    data = {}

    for cell, reads in df:

        segs = df_seg.get_group(cell)

        halfiness_full = []
        halfiness_scaled = []

        for _, row in segs.iterrows():

            seg_chr = row['chr']

            seg_start = row['start']

            seg_end = row['end']

            seg_median = row['median']

            df_seg_bins = reads[(reads['chr'] == seg_chr) &
                                (reads['start'] >= seg_start) &
                                (reads['end'] <= seg_end)]

            df_seg_bins = df_seg_bins[~df_seg_bins['state'].isnull()]

            for state, ideal in zip(df_seg_bins.state, df_seg_bins.ideal):
                if not ideal:
                    continue
                halfiness = - log(np.abs(np.minimum(np.abs(seg_median - state), 0.499) - 0.5), 2) - 1
                scaled = halfiness / (state + 1)

                halfiness_full.append(halfiness)
                halfiness_scaled.append(scaled)

        total_halfiness = np.nansum(halfiness_full)
        scaled_halfiness = np.nansum(halfiness_scaled)

        data[cell] = (total_halfiness, scaled_halfiness)

    return data


def median_of_segment_residuals_from_segment_integer(df_seg):
    '''
    MSRSI: This measures "integerness", but not "dispersion".
    '''
    data = {}

    for cell, segs in df_seg:

        residuals = np.abs(
            segs['median'] -
            segs['state'])

        if len(residuals):
            median_seg_residuals_integer = np.median(residuals)
        else:
            median_seg_residuals_integer = float('nan')

        data[cell] = median_seg_residuals_integer

    return data


def median_of_bin_residuals_from_segment_integer(df):
    '''
    MBRSI: This measures both "dispersion" and "integerness".
    '''
    data = {}

    for cell, reads in df:
        df_hmmcopy = reads[~reads['copy'].isnull()]

        residuals = np.abs(
            df_hmmcopy['copy'] -
            df_hmmcopy['state'])

        median_bin_residuals_integer = np.median(residuals)

        data[cell] = median_bin_residuals_integer
    return data


def compute_total_reads_hmmcopy(df):

    data = {}
    for cell, reads in df:
        total_reads_hmmcopy = sum(reads['reads'][~reads['copy'].isnull()])

        data[cell] = total_reads_hmmcopy

    return data


def compute_mad_hmmcopy(df):

    data = {}
    for cell, reads in df:
        df_hmmcopy = reads[~reads['copy'].isnull()]

        if len(df_hmmcopy) > 0:
            mad_hmmcopy = stand_mad(df_hmmcopy['copy'], c=1)
        else:
            mad_hmmcopy = float('NaN')

        data[cell] = mad_hmmcopy

    return data


def compute_mad_autosomes(df):

    data = {}
    for cell, reads in df:
        df_autosomes = reads[(reads['chr'] != 'X') & (reads['chr'] != 'Y')]

        if len(df_autosomes) > 0:
            mad_autosomes = stand_mad(df_autosomes['copy'], c=1)
        else:
            mad_autosomes = float('NaN')

        data[cell] = mad_autosomes
    return data


def compute_cv_hmmcopy(df):

    data = {}
    for cell, reads in df:

        df_hmmcopy = reads[~reads['copy'].isnull()]

        if len(df_hmmcopy) > 0:
            df_mean = np.mean(df_hmmcopy['copy'])

            df_sd = np.std(df_hmmcopy['copy'])

            cv_hmmcopy = df_sd / df_mean

        else:
            cv_hmmcopy = float('NaN')

        data[cell] = cv_hmmcopy

    return data


def compute_cv_neutral_state(df):

    data = {}
    for cell, reads in df:

        df_neutral_state = reads[
            (reads['state'] == 2) & (
                ~reads['copy'].isnull())]

        if len(df_neutral_state) > 0:
            df_mean = np.mean(df_neutral_state['copy'])

            df_sd = np.std(df_neutral_state['copy'])

            cv_neutral_state = df_sd / df_mean

        else:
            cv_neutral_state = float('NaN')

        data[cell] = cv_neutral_state
    return data


def compute_autocorrelation_hmmcopy(df):

    data = {}
    for cell, reads in df:

        df_hmmcopy = reads[~reads['copy'].isnull()]

        if len(df_hmmcopy) > 0:
            acf_hmmcopy = acf(df_hmmcopy['copy'], nlags=1)

            ac_hmmcopy = acf_hmmcopy[1]

        else:
            ac_hmmcopy = float('NaN')

        data[cell] = ac_hmmcopy
    return data


def compute_mean_median_std_hmmcopy_reads_per_bin(df):

    mean = {}
    median = {}
    std = {}
    for cell, reads in df:

        df_hmmcopy = reads[~reads['copy'].isnull()]

        mean_hmmcopy_reads_per_bin = np.mean(df_hmmcopy['reads'])

        median_hmmcopy_reads_per_bin = np.median(df_hmmcopy['reads'])

        std_hmmcopy_reads_per_bin = np.std(df_hmmcopy['reads'])

        mean[cell] = mean_hmmcopy_reads_per_bin
        median[cell] = median_hmmcopy_reads_per_bin
        std[cell] = std_hmmcopy_reads_per_bin

    return mean, median, std


def compute_num_empty_bins(df):

    data = {}
    for cell, reads in df:

        df_hmmcopy = reads[~reads['copy'].isnull()]

        empty_bins_hmmcopy = len(df_hmmcopy[df_hmmcopy['reads'] == 0])

        data[cell] = empty_bins_hmmcopy

    return data


def compute_mad_neutral_state(df):
    data = {}

    for cell, reads in df:
        df_neutral_state = reads[
            (reads['state'] == 2) & (
                ~reads['copy'].isnull())]

        if len(df_neutral_state) > 0:
            mad_neutral_state = stand_mad(df_neutral_state['copy'], c=1)

        else:
            mad_neutral_state = float('NaN')

        data[cell] = mad_neutral_state

    return data


def median_of_bin_residuals_from_segment_median(df, df_seg):
    '''
    MBRSM: This measures "dispersion", but not "integerness".
    '''

    data = {}

    for cell, reads in df:

        segs = df_seg.get_group(cell)

        residuals = []

        for _, row in segs.iterrows():

            seg_chr = row['chr']

            seg_start = row['start']

            seg_end = row['end']

            seg_median = row['median']

            df_seg_bins = reads[(reads['chr'] == seg_chr) &
                                (reads['start'] >= seg_start) &
                                (reads['end'] <= seg_end)]

            df_seg_bins_hmmcopy = df_seg_bins[~df_seg_bins['copy'].isnull()]

            seg_residuals = [
                np.abs(
                    x -
                    seg_median) for x in df_seg_bins_hmmcopy['copy']]

            residuals.extend(seg_residuals)

        median_bin_residuals_median = np.median(residuals)

        data[cell] = median_bin_residuals_median

    return data
