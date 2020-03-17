import os.path
import sys
from single_cell.tests.jenkins import compare

def compare_count_haps():
    output_path = sys.argv[1]
    ref_path = sys.argv[2]

    refhaps = os.path.join(ref_path, "allele_counts_ref.csv.gz")
    haps = os.path.join(output_path, "allele_counts.csv.gz")

    compare.compare_count_haps(haps, refhaps)

if __name__ == "__main__":
    compare_count_haps()

