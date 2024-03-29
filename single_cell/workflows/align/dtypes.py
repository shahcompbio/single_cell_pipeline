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
        'is_contaminated': 'bool',
        'trim': 'bool',
        'sample_id': 'str',
        'aligned': 'float',
        'expected': 'float',
        'overlap_with_all_filters': 'float',
        'overlap_with_all_filters_and_qual': 'float',
        'overlap_with_dups': 'float',
        'overlap_without_dups': 'float',
        'is_control': 'bool',
    }

    gc = {str(i): 'float' for i in range(0, 101)}
    gc['cell_id'] = 'str'

    dtypes = locals()

    return dtypes


def fastqscreen_dtypes(genome_labels):
    metrics = {'fastqscreen_nohit': 'int', 'cell_id': 'str'}
    for label in genome_labels:
        metrics['fastqscreen_{}'.format(label)] = 'int'
        metrics['fastqscreen_{}_multihit'.format(label)] = 'int'

    fastqscreen_detailed = {
        'cell_id': 'str',
        'readend': 'str',
        'count': 'int'
    }

    for label in genome_labels:
        fastqscreen_detailed[label] = 'int'

    dtypes = locals()
    return dtypes
