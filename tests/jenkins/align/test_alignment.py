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
    metrics += "_alignment_metrics.csv.gz"

    gc_metrics = os.path.join(path, library_id)
    gc_metrics += "_gc_metrics.csv.gz"

    return metrics, gc_metrics

def compare_alignment(ref_metrics, metrics, 
                    ref_gc_metrics, gc_metrics):

    compare.compare_metrics(ref_metrics, metrics)
    compare.compare_metrics(ref_gc_metrics, gc_metrics)

if __name__ == "__main__":

    output_path = sys.argv[1]
    output_lib = sys.argv[2]

    ref_path = sys.argv[3]
    ref_lib = "A97318A"

    ref_metrics, ref_gc_metrics = get_inputs(ref_path, "A97318A")
    metrics, gc_metrics = get_inputs(output_path, output_lib)

    compare_alignment(ref_metrics, metrics,
                    ref_gc_metrics, gc_metrics)
