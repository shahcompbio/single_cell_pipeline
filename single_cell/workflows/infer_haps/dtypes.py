def dtypes():
    haplotypes = {
        'chromosome': 'str',
        'position': 'int',
        'allele': 'str',
        'hap_label': 'str',
        'allele_id': 'str',
        'ref': 'str',
        'alt': 'str'
    }

    dtypes = locals()

    return dtypes
