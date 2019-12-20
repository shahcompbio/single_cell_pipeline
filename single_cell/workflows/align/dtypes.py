def dtypes():
    metrics = {
        'cell_id': 'str',
        'total_mapped_reads': 'Int64',
        'library_id': 'str',
        'unpaired_mapped_reads': 'Int64',
        'paired_mapped_reads': 'Int64',
        'unpaired_duplicate_reads': 'Int64',
        'paired_duplicate_reads': 'Int64',
        'unmapped_reads': 'Int64',
        'percent_duplicate_reads': 'float64',
        'estimated_library_size': 'Int64',
        'total_reads': 'Int64',
        'total_duplicate_reads': 'Int64',
        'total_properly_paired': 'Int64',
        'coverage_breadth': 'float64',
        'coverage_depth': 'float64',
        'median_insert_size': 'float64',
        'mean_insert_size': 'float64',
        'standard_deviation_insert_size': 'float64',
        'cell_call': 'str',
        'column': 'Int64',
        'experimental_condition': 'str',
        'img_col': 'Int64',
        'index_i5': 'str',
        'index_i7': 'str',
        'primer_i5': 'str',
        'primer_i7': 'str',
        'row': 'Int64',
        'sample_type': 'str',
        'fastqscreen_grch37': 'Int64',
        'fastqscreen_salmon': 'Int64',
        'fastqscreen_grch37_multihit': 'Int64',
        'fastqscreen_salmon_multihit': 'Int64',
        'fastqscreen_mm10': 'Int64',
        'fastqscreen_nohit': 'Int64',
        'fastqscreen_mm10_multihit': 'Int64',
        'is_contaminated': 'bool'
    }

    gc = {str(i): 'float64' for i in range(0,100)}
    gc['cell_id'] = 'str'

    dtypes = locals()

    return dtypes
