def dtypes():
    snv_annotate = {
        'cell_id': 'str',
        'chrom': 'str',
        'coord': 'int',
        'ref': 'str',
        'alt': 'str',
        'db_id': 'str',
        'exact_match': 'int',
        'indel': 'int',
        'mappability': 'float',
        'effect': 'str',
        'effect_impact': 'str',
        'functional_class': 'str',
        'codon_change': 'str',
        'amino_acid_change': 'str',
        'amino_acid_length': 'str',
        'gene_name': 'str',
        'transcript_biotype': 'str',
        'gene_coding': 'str',
        'transcript_id': 'str',
        'exon_rank': 'str',
        'genotype': 'str',
        'tri_nucleotide_context': 'str',
    }

    dtypes = locals()

    return dtypes
