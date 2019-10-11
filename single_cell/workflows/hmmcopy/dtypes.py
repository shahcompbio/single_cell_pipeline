def dtypes():
    reads = {
        'chr': 'str',
        'start': 'Int64',
        'end': 'Int64',
        'reads': 'Int64',
        'gc': 'float64',
        'copy': 'float64',
        'state': 'Int64',
        'cell_id': 'str',
        'sample_id': 'str',
        'library_id': 'str',
    }

    segs = {
        'chr': 'str',
        'start': 'Int64',
        'end': 'Int64',
        'state': 'Int64',
        'median': 'float64',
        'multiplier': 'Int64',
        'cell_id': 'str',
        'sample_id': 'str',
        'library_id': 'str',
    }

    metrics = {
        'cell_id': 'str',
        'sample_id': 'str',
        'library_id': 'str',
        'multiplier': 'Int64',
        'MSRSI_non_integerness': 'float64',
        'MBRSI_dispersion_non_integerness': 'float64',
        'MBRSM_dispersion': 'float64',
        'autocorrelation_hmmcopy': 'float64',
        'cv_hmmcopy': 'float64',
        'empty_bins_hmmcopy': 'Int64',
        'mad_hmmcopy': 'float64',
        'mean_hmmcopy_reads_per_bin': 'float64',
        'median_hmmcopy_reads_per_bin': 'Int64',
        'std_hmmcopy_reads_per_bin': 'float64',
        'total_mapped_reads_hmmcopy': 'Int64',
        'total_halfiness': 'float64',
        'scaled_halfiness': 'float64',
        'mean_state_mads': 'float64',
        'mean_state_vars': 'float64',
        'mad_neutral_state': 'float64',
        'breakpoints': 'Int64',
        'mean_copy': 'float64',
        'state_mode': 'Int64',
        'log_likelihood': 'float64',
        'true_multiplier': 'float64',
        'column': 'Int64',
        'img_col': 'Int64',
        'primer_i7': 'str',
        'index_i5': 'str',
        'sample_type': 'str',
        'primer_i5': 'str',
        'experimental_condition': 'str',
        'cell_call': 'str',
        'index_i7': 'str',
        'order': 'Int64',
        'row': 'Int64'
    }

    dtypes = locals()

    return dtypes
