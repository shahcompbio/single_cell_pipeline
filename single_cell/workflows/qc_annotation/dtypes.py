def dtypes():
    metrics = {
        'cell_id': 'category',
        'sample_id': 'category',
        'library_id': 'category',
        'total_mapped_reads_hmmcopy': 'int64',
        'mean_copy': 'float64',
        'state_mode': 'int64',
        'order': 'float64',
        'experimental_condition': 'object',
        'quality': 'float64',
        'is_s_phase': 'bool',
    }

    dtypes = locals()

    return dtypes
