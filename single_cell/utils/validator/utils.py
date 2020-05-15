class DtypeException(Exception):
    pass


class MissingFieldError(Exception):
    pass


class InvalidBarcode(Exception):
    pass


class InvalidIndex(Exception):
    pass


class MissingInput(Exception):
    pass


class InvalidInstrument(Exception):
    pass


class DLPIndexError(Exception):
    pass


def get(data, key):
    if key not in data:
        raise MissingFieldError('{} key missinsg in yaml file.'.format(key))
    return data[key]


def check_data_type(keys, dtype, data):
    for key in keys:

        if not isinstance(get(data, key), dtype):
            raise DtypeException('{} value must be {}'.format(key, dtype))


def check_barcodes(barcode_str):
    for val in barcode_str:
        if val not in ['A', 'C', 'G', 'T']:
            raise InvalidBarcode('{} is not a valid varcode'.format(barcode_str))


def check_sequencing_instrument_type(instrument):
    bam_instruments = ['CAPILLARY', 'DNBSEQ', 'HELICOS', 'ILLUMINA', 'IONTORRENT', 'LS454', 'ONT', 'PACBIO', 'SOLID']

    if instrument not in bam_instruments:
        raise InvalidInstrument(
            '{} instrument is not supported. Please check Bam spec for supported instruments'.format(instrument)
        )



def check_genomic_regions(region, sep='-'):
    chroms = list(map(str,range(1, 23))) + ['X', 'Y']

    chrom, start, end = region.split(sep)

    assert chrom in chroms, '{} is not a valid chrom'.format(chrom)


def check_cells_data(data):
    for cell in data:
        check_data_type(['bam'], str, data[cell])


def check_normal_data(normal):
    if 'bam' in normal:
        check_data_type(['bam'], str, normal)
    else:
        for cell in normal:
            check_data_type(['bam'], str, normal[cell])
