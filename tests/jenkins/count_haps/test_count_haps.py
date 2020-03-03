import os.path
import sys
from single_cell.utils import compare

def compare_count_haps():
    output_path = sys.argv[1]
    ref_path = sys.argv[2]

    refhaps = os.path.join(ref_path, "haplotypes.tsv")
    haps = os.path.join(output_path, "haplotypes.tsv")

    compare.compare_count_haps(haps, refhaps)

if __name__ == "__main__":
    compare_count_haps()

