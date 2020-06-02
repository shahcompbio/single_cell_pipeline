'''
Created on Jun 14, 2018

@author: douglas
'''
import logging
import warnings

import numpy as np
import pandas as pd
from single_cell.utils import csvutils


def _exact_compare_cols(data, reference, column_name):
    data_index = set(data.index)
    reference_index = set(reference.index)

    assert data_index == reference_index

    index_order = sorted(data_index)

    reference = reference.reindex(index_order)
    data = data.reindex(index_order)

    if data[column_name].dtype == float:
        if data[column_name].isnull().all() and reference[column_name].isnull().all():
            return
        assert abs(max(data[column_name] - reference[column_name])) <= 0.001
    else:
        assert data[column_name].equals(reference[column_name])


def _reset_indexes(data, refdata):
    assert data.index.equals(refdata.index)

    index_order = sorted(data.index)
    reindexed_data = data.reindex(index_order)
    reindexed_refdata = refdata.reindex(index_order)

    return reindexed_data, reindexed_refdata


def _approx_compare_cols(data, reference, column_name, eps=0.001):
    data, reference = _reset_indexes(data, reference)

    # check for exact match first
    if data[column_name].equals(reference[column_name]):
        return

    diff = data[column_name] - reference[column_name]

    assert np.nanmax(diff.tolist()) < eps


def _load_hmmcopy_reads_data(readsfile):
    keepcols = [
        'ideal', 'valid', 'gc', 'map', 'state',
        'cor_gc', 'copy'
    ]

    reads = csvutils.read_csv_and_yaml(readsfile)

    reads = reads.set_index(['cell_id', 'chr', 'start', 'end'])

    reads = reads[keepcols]

    return reads


def _load(file, by, reindex=False):
    loaded = csvutils.read_csv_and_yaml(file)

    loaded = loaded.sort_values(by, ascending=[True] * len(by))

    if reindex:
        loaded = loaded.set_index(by)

    return loaded


def _bkps_starts_ends_overlap(data, refdata):
    assert data.chrom1.equals(refdata.chrom1)
    assert data.chrom2.equals(refdata.chrom2)

    for index in refdata.index:
        ref_s1 = refdata.start1[index]
        ref_e1 = refdata.end1[index]
        ref_s2 = refdata.start2[index]
        ref_e2 = refdata.end2[index]

        s1 = data.start1[index]
        e1 = data.end1[index]
        s2 = data.start2[index]
        e2 = data.end2[index]

        if not (ref_s1 <= s1 <= ref_e1) or not (s1 <= ref_e1 <= e1):
            return False

        if not (ref_s2 <= s2 <= ref_e2) or not (s2 <= ref_e2 <= e2):
            return False

    return True


def _call_positions_similar(data, refdata, strict=False,
                            loose=0.005, prob_cutoff=0.55):
    if strict:
        assert refdata.index.equals(data.index)

    shared_calls = refdata.index.intersection(data.index)
    diff_calls = refdata.index.difference(data.index)

    n_shared = shared_calls.size
    percentage_shared = n_shared / refdata.index.size

    if 1 - percentage_shared < loose:
        return True, shared_calls, diff_calls
    return False, shared_calls, diff_calls


def _check_for_missing_cols(data, refdata):
    data_cols = set(data.columns.values)
    ref_cols = set(refdata.columns.values)

    common_cols = list(data_cols.intersection(ref_cols))

    missing_cols = list(ref_cols - data_cols)

    if missing_cols:
        warnings.warn("missing cols in a reference: {}".format(missing_cols))

    return common_cols


def compare_count_haps(haps, refhaps):
    haps = pd.read_csv(haps)
    refhaps = pd.read_csv(refhaps)

    if haps.empty and refhaps.empty:
        logging.getLogger('testing').warning("comparing empty hap counts")
        return

    cols_must_match = ["chromosome", "start", "end", "allele_id", "hap_label", "cell_id", "readcount"]

    for col in cols_must_match:
        _exact_compare_cols(haps, refhaps, col)


def compare_infer_haps(data, refdata):
    data = pd.read_csv(data)
    refdata = pd.read_csv(refdata)

    assert data.equals(refdata)


def compare_tables(data, refdata, eps=None):
    common_cols = _check_for_missing_cols(data, refdata)
    for colname in common_cols:
        if not eps:
            _exact_compare_cols(data, refdata, colname)
        else:
            _approx_compare_cols(data, refdata, colname, eps=eps)


def compare_annotation(annotation, refannotation):
    annotation = csvutils.read_csv_and_yaml(annotation)
    refannotation = csvutils.read_csv_and_yaml(refannotation)

    common_cols = _check_for_missing_cols(annotation, refannotation)
    for col in common_cols:
        ann = annotation[col].dropna()
        ref = refannotation[col].dropna()
        assert set(ann) == set(ref)


def compare_variant_calls(callsdata, refcallsdata):
    by = ["chrom", "coord"]

    calls = _load(callsdata, by, reindex=True)
    refcalls = _load(refcallsdata, by, reindex=True)

    assert calls.index.size > 0
    assert refcalls.index.size > 0

    similar, shared, diff = _call_positions_similar(calls, refcalls)

    assert similar

    diff_data = calls[calls.index.isin(diff)]

    if diff_data.empty:
        return

    assert (diff_data['score'] < 0.55).all()

    calls = calls[calls.index.isin(shared)]
    refcalls = refcalls[refcalls.index.isin(shared)]

    compare_tables(calls, refcalls)


def compare_breakpoint_calls(calls, refcalls):
    sort_by = ["breakpoint_id"]

    calls = _load(calls, sort_by, True)
    refcalls = _load(refcalls, sort_by, True)

    assert _bkps_starts_ends_overlap(calls, refcalls)

    exact_match = ["strands", "type"]
    approx_match = ["max_pos1", "max_chr1", "max_pos1", "confidence_interval_chr1",
                    "confidence_interval_start1", "confidence_interval_end1", "max_chr2",
                    "max_pos2", "confidence_interval_chr2", "confidence_interval_start2",
                    "confidence_interval_end2", "score", "tumour_SR", "tumour_PE",
                    "normal_SR", "normal_PE"]

    assert all(col in refcalls.columns for col in calls.columns)

    for col in exact_match:
        _exact_compare_cols(calls, refcalls, col)
    for col in approx_match:
        _approx_compare_cols(calls, refcalls, col)


def compare_reads(readsdata, refreadsdata):
    reads = _load_hmmcopy_reads_data(readsdata)
    refreads = _load_hmmcopy_reads_data(refreadsdata)

    _exact_compare_cols(reads, refreads, 'ideal')
    _exact_compare_cols(reads, refreads, 'valid')
    _exact_compare_cols(reads, refreads, 'state')
    _exact_compare_cols(reads, refreads, 'gc')
    _exact_compare_cols(reads, refreads, 'map')
    _approx_compare_cols(reads, refreads, 'copy')


def compare_metrics(metrics, refmetrics):
    metrics = _load(metrics, ["cell_id"], reindex=True)
    refmetrics = _load(refmetrics, ["cell_id"], reindex=True)

    compare_tables(metrics, refmetrics)


def compare_annotation_metrics(metrics, refmetrics):
    metrics = _load(metrics, ['cell_id'], reindex=True)
    refmetrics = _load(refmetrics, ['cell_id'], reindex=True)

    exact_cols = [
        'multiplier', 'MSRSI_non_integerness', 'MBRSI_dispersion_non_integerness', 'empty_bins_hmmcopy',
        'total_mapped_reads_hmmcopy', 'breakpoints', 'state_mode', 'column', 'img_col', 'row', 'order',
        'fastqscreen_mm10', 'unpaired_duplicate_reads', 'paired_mapped_reads', 'estimated_library_size',
        'paired_duplicate_reads', 'fastqscreen_salmon_multihit', 'total_properly_paired',
        'total_duplicate_reads', 'fastqscreen_salmon', 'unpaired_mapped_reads', 'unmapped_reads',
        'fastqscreen_grch37_multihit', 'total_reads', 'fastqscreen_nohit', 'total_mapped_reads',
        'fastqscreen_mm10_multihit', 'fastqscreen_grch37',
    ]
    compare_tables(metrics[exact_cols], refmetrics[exact_cols])

    approx_cols = [
        'MSRSI_non_integerness', 'MBRSI_dispersion_non_integerness', 'MBRSM_dispersion',
        'autocorrelation_hmmcopy', 'cv_hmmcopy', 'mad_hmmcopy', 'mean_hmmcopy_reads_per_bin',
        'median_hmmcopy_reads_per_bin', 'std_hmmcopy_reads_per_bin', 'total_halfiness', 'scaled_halfiness',
        'mean_state_mads', 'mean_state_vars', 'mad_neutral_state', 'mean_copy', 'log_likelihood',
        'true_multiplier', 'is_s_phase_prob', 'quality', 'percent_duplicate_reads', 'mean_insert_size',
        'coverage_breadth', 'standard_deviation_insert_size', 'coverage_depth', 'median_insert_size',
    ]

    compare_tables(metrics[approx_cols], refmetrics[approx_cols], eps=0.02)
