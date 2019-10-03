
def dtypes():
    reads = {
        'chr': 'category',
        'start': 'int64',
        'end': 'int64',
        'reads': 'int64',
        'gc': 'float64',
        'copy': 'float64',
        'state': 'int64',
        'cell_id': 'category',
        'sample_id': 'category',
        'library_id': 'category',
    }

    segs = {
        'chr': 'category',
        'start': 'int64',
        'end': 'int64',
        'state': 'int64',
        'median': 'float64',
        'multiplier': 'int64',
        'cell_id': 'category',
        'sample_id': 'category',
        'library_id': 'category',
    }

    metrics = {
        'cell_id': 'category',
        'sample_id': 'category',
        'library_id': 'category',
        'multiplier': 'int64',
        'MSRSI_non_integerness': 'float64',
        'MBRSI_dispersion_non_integerness': 'float64',
        'MBRSM_dispersion': 'float64',
        'autocorrelation_hmmcopy': 'float64',
        'cv_hmmcopy': 'float64',
        'empty_bins_hmmcopy': 'int64',
        'mad_hmmcopy': 'float64',
        'mean_hmmcopy_reads_per_bin': 'float64',
        'median_hmmcopy_reads_per_bin': 'int64',
        'std_hmmcopy_reads_per_bin': 'float64',
        'total_mapped_reads_hmmcopy': 'int64',
        'total_halfiness': 'float64',
        'scaled_halfiness': 'float64',
        'mean_state_mads': 'float64',
        'mean_state_vars': 'float64',
        'mad_neutral_state': 'float64',
        'breakpoints': 'int64',
        'mean_copy': 'float64',
        'state_mode': 'int64',
        'log_likelihood': 'float64',
        'true_multiplier': 'float64',
        'column': 'int64',
        'img_col': 'int64',
        'primer_i7': 'str',
        'index_i5': 'str',
        'sample_type': 'str',
        'primer_i5': 'str',
        'experimental_condition': 'category',
        'cell_call': 'str',
        'index_i7': 'str',
        'order': 'int64',
        'row': 'int64'        
    }

    dtypes = locals()

    return dtypes
