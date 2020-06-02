import os.path
import sys

from single_cell.tests.jenkins import compare


def compare_infer_haps():
    output_path = sys.argv[1]
    ref_path = sys.argv[2]

    refhaps = os.path.join(ref_path, "ref_haplotypes.csv.gz")
    haps = os.path.join(output_path, "haplotypes.csv.gz")

    compare.compare_infer_haps(haps, refhaps)


if __name__ == "__main__":
    compare_infer_haps()
