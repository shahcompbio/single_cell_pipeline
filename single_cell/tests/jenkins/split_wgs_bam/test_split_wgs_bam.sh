import os.path
import sys
from single_cell.utils import compare


def compare_split_counts():
    output_path = sys.argv[1]
    ref_path = sys.argv[2]

    refcounts = os.path.join(ref_path, "counts.csv")
    counts = os.path.join(counts, "counts.csv")

  compare.compare_metrics(refcounts, counts)

if __name__ == "__main__":
    compare_split_counts()

