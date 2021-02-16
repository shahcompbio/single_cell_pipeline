def dtypes():
    metrics = {
        'cell_id': 'str',
        'total_mapped_reads': 'int',
        'library_id': 'str',
        'unpaired_mapped_reads': 'int',
        'paired_mapped_reads': 'int',
        'unpaired_duplicate_reads': 'int',
        'paired_duplicate_reads': 'int',
        'unmapped_reads': 'int',
        'percent_duplicate_reads': 'float',
        'estimated_library_size': 'int',
        'total_reads': 'int',
        'total_duplicate_reads': 'int',
        'total_properly_paired': 'int',
        'coverage_breadth': 'float',
        'coverage_depth': 'float',
        'median_insert_size': 'float',
        'mean_insert_size': 'float',
        'standard_deviation_insert_size': 'float',
        'cell_call': 'str',
        'column': 'int',
        'experimental_condition': 'str',
        'img_col': 'int',
        'index_i5': 'str',
        'index_i7': 'str',
        'primer_i5': 'str',
        'primer_i7': 'str',
        'row': 'int',
        'sample_type': 'str',
        'fastqscreen_grch37': 'int',
        'fastqscreen_salmon': 'int',
        'fastqscreen_grch37_multihit': 'int',
        'fastqscreen_salmon_multihit': 'int',
        'fastqscreen_mm10': 'int',
        'fastqscreen_nohit': 'int',
        'fastqscreen_mm10_multihit': 'int',
        'is_contaminated': 'bool',
        'trim': 'bool',
        'library_id': 'str',
        'sample_id':'str'
    }

    gc = {str(i): 'float' for i in range(0,101)}
    gc['cell_id'] = 'str'

    fastqscreen_detailed = {
        'cell_id': 'str',
        'readend': 'str',
        'grch37': 'int',
        'mm10': 'int',
        'salmon': 'int',
        'count': 'int'
    }

    dtypes = locals()

    return dtypes
