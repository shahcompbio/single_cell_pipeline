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
        'order': 'float64',
        'total_mapped_reads_hmmcopy': 'int64',
        'experimental_condition': 'object',
        'mean_copy': 'float64',
        'state_mode': 'int64',
    }

    dtypes = locals()

    return dtypes
