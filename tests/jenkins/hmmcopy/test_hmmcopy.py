import os.path
import sys
from single_cell.utils import compare


def get_inputs(path, library_id):
    '''
    get metrics and gc metrics given a directory and library
    :param path:  path to metrics files
    :param library_id: library id associated with metrics files
    '''
    metrics = os.path.join(path, library_id)
    metrics += "_hmmcopy_metrics.csv.gz"

    reads = os.path.join(path, library_id)
    reads += "_reads.csv.gz"

    return metrics, reads


if __name__ == "__main__":
    output_path = sys.argv[1]
    output_lib = sys.argv[2]

    ref_path = sys.argv[3]
    ref_lib = "A97318A"

    ref_metrics, ref_reads = get_inputs(ref_path, "A97318A")
    metrics, reads = get_inputs(output_path, output_lib)

    compare.compare_metrics(ref_metrics, metrics)
    compare.compare_reads(ref_reads, reads)