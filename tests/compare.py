'''
Created on Jun 14, 2018

@author: dgrewal
'''
import pandas as pd
import numpy as np
import argparse


def exact_compare_cols(data, reference, column_name):

    assert data.index.equals(reference.index)

    assert data[column_name].equals(reference[column_name])


def approx_compare_cols(data, reference, column_name, eps=0.001):

    assert data.index.equals(reference.index)

    diff = data[column_name] - reference[column_name]

    assert np.nanmax(diff) < eps


def read_hdf(tablename, filename, keep_cols=None):

    df = pd.read_hdf(filename, key=tablename, usecols=keep_cols)

    return df


def read_hdf_chunks(tablename, filename, keep_cols=None):
    data = []
    chunksize = 10 ** 6

    for chunk in pd.read_hdf(
            filename, chunksize=chunksize, key=tablename, usecols=keep_cols):
        data.append(chunk)

    df = pd.concat(data)

    return df


def load_hmmcopy_reads_data(readsfile, tablename):
    keepcols = [
        'ideal', 'valid', 'gc', 'map', 'state',
        'cor_gc', 'copy']

    reads = read_hdf_chunks(tablename, readsfile, keep_cols=keepcols)

    reads = reads.set_index(['cell_id', 'chr', 'start', 'end'])

    return reads


def load_metrics_data(filename, tablename):

    reads = read_hdf(tablename, filename, keep_cols=None)

    reads = reads.set_index(['cell_id'])

    return reads


def compare_tables(data, refdata):

    assert np.array_equal(data.columns.values, refdata.columns.values)

    for colname in data.columns.values:
        exact_compare_cols(data, refdata, colname)


def compare_hmmcopy(readsdata, refreadsdata):

    for multiplier in range(7):
        tablename = '/hmmcopy/reads/{}'.format(multiplier)

        reads = load_hmmcopy_reads_data(readsdata, tablename)
        refreads = load_hmmcopy_reads_data(refreadsdata, tablename)

        exact_compare_cols(reads, refreads, 'ideal')
        exact_compare_cols(reads, refreads, 'valid')
        exact_compare_cols(reads, refreads, 'state')
        exact_compare_cols(reads, refreads, 'gc')
        exact_compare_cols(reads, refreads, 'map')

        approx_compare_cols(reads, refreads, 'copy')

        tablename = '/hmmcopy/metrics/{}'.format(multiplier)
        metrics = load_metrics_data(readsdata, tablename)
        refmetrics = load_metrics_data(refreadsdata, tablename)

        compare_tables(metrics, refmetrics)


def compare_alignment(metricsdata, refmetricsdata):

    tablename = '/alignment/metrics/'

    metrics = load_metrics_data(metricsdata, tablename)
    refmetrics = load_metrics_data(refmetricsdata, tablename)

    compare_tables(metrics, refmetrics)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--hmmcopy')
    parser.add_argument('--alignment')

    parser.add_argument('--ref_hmmcopy')
    parser.add_argument('--ref_alignment')

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = parse_args()

    compare_alignment(args.alignment, args.ref_alignment)

    compare_hmmcopy(args.hmmcopy, args.ref_hmmcopy)

