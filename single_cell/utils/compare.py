'''
Created on Jun 14, 2018

@author: dgrewal
'''
import warnings
import numpy as np
from single_cell.utils import csvutils
import pysam
import os
import pandas as pd

def exact_compare_cols(data, reference, column_name):
    data_index = set(data.index)
    reference_index = set(reference.index)

    assert data_index == reference_index

    index_order = sorted(data_index)
    
    reference = reference.reindex(index_order)
    data = data.reindex(index_order)

    assert data[column_name].equals(reference[column_name])

def reset_indexes(data, refdata):

    assert data.index.equals(refdata.index)

    index_order = sorted(data.index)
    reindexed_data = data.reindex(index_order)
    reindexed_refdata = refdata.reindex(index_order)

    return reindexed_data, reindexed_refdata
    
def approx_compare_cols(data, reference, column_name, eps=0.001):
    data_index = set(data.index)
    reference_index = set(reference.index)

    data, reference = reset_indexes(data, reference)

    #check for exact match first
    if data[column_name].equals(reference[column_name]):
        return

    diff = data[column_name] - reference[column_name]

    assert np.nanmax(diff.tolist()) < eps

def load_hmmcopy_reads_data(readsfile):
    keepcols = [
        'ideal', 'valid', 'gc', 'map', 'state',
        'cor_gc', 'copy']

    reads = csvutils.read_csv_and_yaml(readsfile)

    reads = reads.set_index(['cell_id', 'chr', 'start', 'end'])

    reads = reads[keepcols]

    return reads

def load(file, by, reindex = False):

    loaded = csvutils.read_csv_and_yaml(file)

    loaded = loaded.sort_values(by, ascending = [True] * len(by))

    if reindex:
        loaded = loaded.set_index(by)

    return loaded

def _all_low_PR(calls, probname):
    is_low_PR = calls[probname] < 0.55
    is_low_PR = is_low_PR.tolist()

    return not False in is_low_PR

def has_match(data, s1, e1, s2, e2):
    matches_s1 = data.loc[(data.start1 >= s1 ) & (data.start1 <= e1 )]
    matches_all = data.loc[(data.start2 >= s2 ) & (data.start2 <= e2 )]

    if matches_all.index.size >= 1: 
        return True
    return False

def bkps_starts_ends_overlap(data, refdata):
    print (data)
    assert set(data.chrom1.unique()) == set(refdata.chrom1.unique())
    assert set(data.chrom2.unique()) == set(refdata.chrom2.unique())

    chromsets = set([(c1, c2) for c1, c2 in zip(data.chrom1, data.chrom2)])

    for chroms in chromsets:
        chromdata = data.loc[(data.chrom1 == chroms[0]) 
                                & (data.chrom2 == chroms[1])]
        chromref = refdata.loc[(refdata.chrom1 == chroms[0]) 
                                & (refdata.chrom2 == chroms[1])]

        for (i, refbkt), (i, bkt) in zip(chromref.iterrows(), chromdata.iterrows()):
            if not has_match(chromdata, refbkt.start1, refbkt.end1, 
                                refbkt.start2, refbkt.end2):
                return False

            if not has_match(chromref, bkt.start1, bkt.end1, 
                                bkt.start2, bkt.end2):
                return False
    return True
        
def compare_infer_haps(data, refdata):
    data = load(data, ["chromosome", "position"], reindex = True)
    refdata = load(refdata, ["chromosome", "position"], reindex = True)

    similar, shared, diff = call_positions_similar(data, refdata)

    assert similar

    data = data[data.index.isin(shared_calls)]
    refdata = refdata[refdata.index.isin(shared_calls)]

    compare_tables(data, refdata)

def call_positions_similar(data, refdata, strict = False, 
                        loose = 0.005, prob_cutoff = 0.55):
    if strict:
        assert refdata.index.equals(data.index)

    shared_calls = refdata.index.intersection(data.index)
    diff_calls = refdata.index.difference(data.index)

    n_shared = shared_calls.size
    percentage_shared = n_shared/refdata.index.size

    if 1 - percentage_shared < loose:
        return True, shared_calls, diff_calls
    return False, shared_calls, diff_calls

def check_for_missing_cols(data, refdata):
    data_cols = set(data.columns.values)
    ref_cols = set(refdata.columns.values)

    common_cols = list(data_cols.intersection(ref_cols))

    missing_cols = list(ref_cols - data_cols)

    if missing_cols:
        warnings.warn("missing cols in a reference: {}".format(missing_cols))
    
    return common_cols

def compare_tables(data, refdata):
    common_cols = check_for_missing_cols(data, refdata)
    for colname in common_cols:
        exact_compare_cols(data, refdata, colname)

def compare_annotation(annotation, refannotation):
    annotation = csvutils.read_csv_and_yaml(annotation)
    refannotation = csvutils.read_csv_and_yaml(refannotation)

    common_cols = check_for_missing_cols(annotation, refannotation)
    for col in common_cols:
        ann = annotation[col].dropna()
        ref = refannotation[col].dropna()
        assert set(ann) == set(ref)

def compare_variant_calls(callsdata, refcallsdata):
    by = ["chrom", "coord"]

    calls = load(callsdata, by, reindex = True)
    refcalls = load(refcallsdata, by, reindex = True)

    assert calls.index.size > 0 
    assert refcalls.index.size > 0

    similar, shared, diff = call_positions_similar(calls, refcalls)

    assert similar

    diff_data = calls[calls.index.isin(diff)] 

    assert _all_low_PR(diff_data, "score")

    calls = calls[calls.index.isin(shared)]
    refcalls = refcalls[refcalls.index.isin(shared)]

    compare_tables(calls, refcalls)

def overlap(start1, end1, start2, end2):
    """Does the range (start1, end1) overlap with (start2, end2)?"""
    return end1 >= start2 and end2 >= start1

def readcount_matches_target(mapped_target, unmapped_target, 
                            bam, loose = 0.0001):
    mapped = reduce(lambda x, y: x + y, [int(l.rstrip('\n').split('\t')[2]) for l in pysam.idxstats(bam)])
    
    unmapped = reduce(lambda x, y: x + y, 
        [int(l.rstrip('\n').split('\t')[3]) for l in pysam.idxstats(bam)])
    
    allowance_mapped = mapped * loose
    allowance_unmapped = unmapped * loose

    mapped_range = [mapped - allowance_mapped, mapped + allowance_mapped]
    unmapped_range = [unmapped - allowance_unmapped, unmapped + allowance_unmapped]

    mapped_matches = mapped_target >= mapped_range[0] 
                    and mapped_target <= mapped_range[1] 

    unmapped_matches = unmapped_target >= unmapped_range[0] 
                    and unmapped_target <= unmapped_range[1] 
    
        
    return mapped_matches and unmapped_matches

def parse_bams_from_dir(dir):
    
    onlyfiles = [f for f in os.listdir(dir) 
                    if os.path.isfile(os.path.join(dir, f))]

    onlyfiles = [f for f in onlyfiles if not ".bai" in f ]

    onlyfiles = [f for f in onlyfiles if  ".bam" in f ]

    intervals = [get_interval_from_bamname(f) for f in onlyfiles]

    onlyfiles = [os.path.join(dir, f) for f in onlyfiles]

    return intervals, onlyfiles

def get_interval_from_bamname(bam):
    return bam.split(".")[0]

def assess_bam_files(compare_dir, expected_counts):

    assert os.path.isdir(compare_dir)
    
    expected_counts = pd.read_csv(expected_counts, sep = ",")
    intervals, bams = parse_bams_from_dir(compare_dir)

    for interval, bam in zip(intervals, bams):
        assert bam_counts_match(expected_counts, interval, bam)
        assert bam_header_wellformed(bam)

def bam_header_wellformed(bam):
    bam = pysam.AlignmentFile(bam, "rb")
    header = bam.text  

    if not header: #fails when empty
        return False 
    if not "@PG" in header:
        return False
    if not "@RG" in header:
        return False
    if not "@SQ" in header:
        return False  
        
    return True

def bam_counts_match(expected_counts, interval, bam):

    at_interval = expected_counts.loc[expected_counts.interval == interval]
    at_interval = at_interval.to_dict(orient="list")

    target_mapped = at_interval["mapped"][0]
    target_unmapped = at_interval["unmapped"][0]

    if not readcount_matches_target(target_mapped, target_unmapped, bam):
        return False
    return True

def compare_breakpoint_calls(calls, refcalls): 

    sort_by = ["chrom1", "start1", "chrom2", "start2"]

    calls = load(calls, sort_by)
    refcalls = load(refcalls, sort_by)

    assert bkps_starts_ends_overlap(calls, refcalls)

    exact_match = ["strands", "type"]
    approx_match = ["max_pos1", "max_chr1", "max_pos1", "confidence_interval_chr1",
                    "confidence_interval_start1", "confidence_interval_end1", "max_chr2",
                    "max_pos2", "confidence_interval_chr2","confidence_interval_start2", 
                    "confidence_interval_end2", "score", "None_SR", "None_PE"]

    common_cols = check_for_missing_cols(calls, refcalls)

    assert all(col in common_cols for col in approx_match + exact_match)

    for col in exact_match:
        exact_compare_cols(calls, refcalls, col)
    for col in approx_match:
        approx_compare_cols(calls, refcalls, col)


def compare_reads(readsdata, refreadsdata):
    reads = load_hmmcopy_reads_data(readsdata)
    refreads = load_hmmcopy_reads_data(refreadsdata)

    exact_compare_cols(reads, refreads, 'ideal')
    exact_compare_cols(reads, refreads, 'valid')
    exact_compare_cols(reads, refreads, 'state')
    exact_compare_cols(reads, refreads, 'gc')
    exact_compare_cols(reads, refreads, 'map')
    approx_compare_cols(reads, refreads, 'copy')


def compare_metrics(metrics, refmetrics):
    metrics = load(metrics, ["cell_id"], reindex = True)
    refmetrics = load(refmetrics, ["cell_id"], reindex = True)

    compare_tables(metrics, refmetrics)

