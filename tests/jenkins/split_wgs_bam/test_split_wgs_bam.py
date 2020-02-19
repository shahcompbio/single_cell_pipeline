import os.path
import sys
from single_cell.utils import compare
import pandas as pd

def compare_split_counts():
    output_path = sys.argv[1]
    ref_path = sys.argv[2]

    refcounts = os.path.join(ref_path, "counts.csv")
    counts = os.path.join(output_path, "counts.csv")

    counts = pd.read_csv(counts)
    refcounts = pd.read_csv(refcounts)

    compare.compare_tables(counts, refcounts)


if __name__ == "__main__":
    compare_split_counts()

