import os

from single_cell.utils import inpututils


def create_genotyping_workflow(args):
    strelka_vcf, museq_vcf, tumour_cell_bams = inpututils.load_variant_counting_input(
        args['input_yaml']
    )

    counts_output = os.path.join(args['out_dir'], "counts.csv.gz")

    config = inpututils.load_config(args)
    config = config['variant_calling']


def genotyping_pipeline(args):
    raise NotImplementedError()
