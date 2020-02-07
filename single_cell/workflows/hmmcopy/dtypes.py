def dtypes():
    reads = {
        'chr': 'str',
        'start': 'int',
        'end': 'int',
        'width': 'int',
        'reads': 'int',
        'gc': 'float',
        'cor_gc': 'float',
        'cor_map': 'float',
        'copy': 'float',
        'map': 'float',
        'state': 'int',
        'cell_id': 'str',
        'sample_id': 'str',
        'library_id': 'str',
        'valid': 'bool',
        'ideal': 'bool',
        'modal_curve': 'float',
        'modal_quantile': 'float',
        'multiplier': 'int',
        'is_low_mappability': 'bool'
    }

    segs = {
        'chr': 'str',
        'start': 'int',
        'end': 'int',
        'state': 'int',
        'median': 'float',
        'multiplier': 'int',
        'cell_id': 'str',
    }

    params = {
        'iteration': 'float',
        # 'is_final': 'bool',
        'state':'float',
        'parameter': 'str',
        'cell_id':'str',
        'value':'float',
    }

    metrics = {
        'multiplier': 'int',
        'cell_id': 'str',
        'sample_id': 'str',
        'library_id': 'str',
        'MSRSI_non_integerness': 'float',
        'MBRSI_dispersion_non_integerness': 'float',
        'MBRSM_dispersion': 'float',
        'autocorrelation_hmmcopy': 'float',
        'cv_hmmcopy': 'float',
        'empty_bins_hmmcopy': 'int',
        'mad_hmmcopy': 'float',
        'mean_hmmcopy_reads_per_bin': 'float',
        'median_hmmcopy_reads_per_bin': 'float',
        'std_hmmcopy_reads_per_bin': 'float',
        'total_mapped_reads_hmmcopy': 'int',
        'total_halfiness': 'float',
        'scaled_halfiness': 'float',
        'mean_state_mads': 'float',
        'mean_state_vars': 'float',
        'mad_neutral_state': 'float',
        'breakpoints': 'int',
        'mean_copy': 'float',
        'state_mode': 'int',
        'log_likelihood': 'float',
        'true_multiplier': 'float',
        'column': 'int',
        'img_col': 'int',
        'primer_i7': 'str',
        'index_i5': 'str',
        'sample_type': 'str',
        'primer_i5': 'str',
        'experimental_condition': 'str',
        'cell_call': 'str',
        'index_i7': 'str',
        'order': 'int',
        'row': 'int'
    }

    dtypes = locals()

    return dtypes
