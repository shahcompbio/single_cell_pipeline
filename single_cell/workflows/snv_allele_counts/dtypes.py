def dtypes():
    snv_allele_counts = {
        'chrom': 'str',
        'coord': 'int',
        'ref': 'str',
        'alt': 'str',
        'ref_counts': 'int',
        'alt_counts': 'int',
        'cell_id': 'str',
        'sample_id': 'str',
        'library_id': 'str',
    }

    dtypes = locals()

    return dtypes
