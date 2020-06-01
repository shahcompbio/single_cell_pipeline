def dtypes():
    readcount = {
        'chromosome': 'str',
        'start': 'int',
        'end': 'int',
        'hap_label': 'str',
        'allele_id': 'str',
        'readcount': 'int',
        'cell_id': 'str'
    }

    dtypes = locals()

    return dtypes
