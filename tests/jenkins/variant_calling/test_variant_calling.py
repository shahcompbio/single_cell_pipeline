import os.path
import sys
from single_cell.utils import compare
from single_cell.utils import csvutils

def get_inputs(path):
    """"
    get metrics and gc metrics given a directory and library
    :param path:  path to metrics files
    """
    strelka = os.path.join(path, "snv_strelka.csv.gz")
    museq = os.path.join(path, "snv_museq.csv.gz")
    snpeff = os.path.join(path, "snv_snpeff.csv.gz")

    return strelka, museq, snpeff


def test_breakpoint_calling(args):
    output_path = args[1]
    ref_path = args[2]

    ref_strelka, ref_museq, ref_snpeff = get_inputs(ref_path)
    strelka, museq, snpeff = get_inputs(output_path)

    # compare.compare_variant_calls(ref_museq, museq)
    compare.compare_variant_calls(ref_snpeff, snpeff)

    ref_strelka = csvutils.read_csv_and_yaml(ref_strelka)
    strelka = csvutils.read_csv_and_yaml(strelka)

    assert ref_strelka.empty and strelka.empty

if __name__ == "__main__":
    test_breakpoint_calling(sys.argv)
