from collections import OrderedDict

import pandas as pd
import vcf
from bx.bbi.bigwig_file import BigWigFile
from single_cell.utils import csvutils
from single_cell.workflows.mappability_annotation.dtypes import dtypes


def parse_region_for_vcf(region):
    if ':' not in region:
        return region, None, None

    chrom, coords = region.split(':')

    if '-' not in coords:
        return chrom, int(coords) - 1, None

    beg, end = coords.split('-')

    beg = int(beg) - 1
    end = int(end)

    return chrom, beg, end


def get_mappability(
        mappability_file,
        vcf_file,
        out_file,
        region=None,
        append_chr=True):
    map_reader = BigWigFile(open(mappability_file, 'rb'))

    vcf_reader = vcf.Reader(filename=vcf_file)

    if region is not None:
        chrom, beg, end = parse_region_for_vcf(region)
        try:
            vcf_reader = vcf_reader.fetch(chrom, start=beg, end=end)
        except ValueError:
            print("no data for region {} in vcf".format(region))
            vcf_reader = []

    data = []

    for record in vcf_reader:
        if append_chr:
            chrom = 'chr{0}'.format(record.CHROM)

        else:
            chrom = record.CHROM

        coord = record.POS

        beg = coord - 100

        beg = max(beg, 0)

        end = coord + 100

        result = map_reader.query(chrom, beg, end, 1)

        if result is None:
            mappability = 0

        else:
            mappability = result[0]['mean']

        data.append({'chrom': record.CHROM, 'coord': record.POS, 'mappability': mappability})

    data = pd.DataFrame(data)

    csvutils.write_dataframe_to_csv_and_yaml(data, out_file, dtypes())


def get_regions(chromosome_lengths, split_size):
    if split_size is None:
        return dict(enumerate(chromosome_lengths.keys()))

    regions = {}
    region_index = 0

    for chrom, length in chromosome_lengths.items():
        lside_interval = range(1, length + 1, split_size)
        rside_interval = range(split_size, length + split_size, split_size)

        for beg, end in zip(lside_interval, rside_interval):
            end = min(end, length)

            regions[region_index] = '{}:{}-{}'.format(chrom, beg, end)
            region_index += 1

    return regions


def get_vcf_regions(vcf_file, split_size):
    chromosome_lengths = load_vcf_chromosome_lengths(vcf_file)
    return get_regions(chromosome_lengths, split_size)


def calculate_vcf_chromosome_lengths(file_name):
    if file_name.endswith('gz'):
        compression = 'gzip'
    else:
        compression = None

    tsv_reader = pd.read_csv(
        file_name, sep='\t', comment='#', chunksize=int(1e6),
        names=['chrom', 'coord'], usecols=[0, 1], index_col=0,
        converters={'chrom': str}, compression=compression)

    max_coord = [chunk.groupby(level=0).max() for chunk in tsv_reader]
    max_coord = pd.concat(max_coord).groupby(level=0).max()

    chromosome_lengths = max_coord + 1000

    return chromosome_lengths['coord'].to_dict()


def load_vcf_chromosome_lengths(file_name, chromosomes=None):
    chromosome_lengths = OrderedDict()

    vcf_reader = vcf.Reader(filename=file_name)

    if len(vcf_reader.contigs) == 0:
        return calculate_vcf_chromosome_lengths(file_name)

    if chromosomes is None:
        chromosomes = vcf_reader.contigs.keys()

    else:
        chromosomes = chromosomes

    calc_lens = calculate_vcf_chromosome_lengths(file_name)

    for chrom, contig in vcf_reader.contigs.items():
        assert chrom == contig.id

        if chrom not in chromosomes:
            continue

        if contig.length is None:
            chromosome_lengths[str(chrom)] = calc_lens[str(chrom)]

        else:
            chromosome_lengths[str(chrom)] = int(contig.length)

    if len(chromosome_lengths) == 0:
        raise Exception('no chromosomes found in vcf header')

    return chromosome_lengths
