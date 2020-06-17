def dtypes():
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
        'median_hmmcopy_reads_per_bin': 'float64',
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
        'row': 'Int64',
        'is_s_phase': 'bool',
        'is_s_phase_prob': 'float64',
        'quality': 'float64',
        'fastqscreen_mm10': 'Int64',
        'coverage_depth': 'float64',
        'paired_duplicate_reads': 'Int64',
        'fastqscreen_grch37': 'Int64',
        'total_reads': 'Int64',
        'unpaired_duplicate_reads': 'Int64',
        'percent_duplicate_reads': 'float64',
        'coverage_breadth': 'float64',
        'mean_insert_size': 'float64',
        'unpaired_mapped_reads': 'Int64',
        'fastqscreen_grch37_multihit': 'Int64',
        'median_insert_size': 'float64',
        'fastqscreen_salmon': 'Int64',
        'total_duplicate_reads': 'Int64',
        'is_contaminated': 'bool',
        'fastqscreen_mm10_multihit': 'Int64',
        'fastqscreen_salmon_multihit': 'Int64',
        'estimated_library_size': 'Int64',
        'standard_deviation_insert_size': 'float64',
        'fastqscreen_nohit': 'Int64',
        'unmapped_reads': 'Int64',
        'total_mapped_reads': 'Int64',
        'total_properly_paired': 'Int64',
        'paired_mapped_reads': 'Int64',
        'order_corrupt_tree': 'Int64',
        'species': 'str',
    }

    dtypes = locals()

    return dtypes
