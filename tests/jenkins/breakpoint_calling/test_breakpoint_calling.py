import os.path 
import sys 
from single_cell.utils import compare

def get_inputs(path, lib = ""):
    '''
    get metrics and gc metrics given a directory and library
    :param path:  path to metrics files
    :param library_id: library id associated with metrics files
    '''

    if lib:
        lib = "_" + lib
    
    lumpy = os.path.join(path, lib + "lumpy_breakpoints.csv.gz") 

    destruct = os.path.join(path, lib + "destruct_breakpoints.csv.gz") 

    return lumpy, destruct 

def compare_breakpoint_calling(ref_lumpy, lumpy, 
                            ref_destruct, destruct):

        compare.compare_breakpoint_calls(ref_lumpy, lumpy)
        compare.compare_breakpoint_calls(ref_destruct, destruct)

if __name__ == "__main__":
    output_path = sys.argv[1]
    ref_path = sys.argv[2]

    ref_lumpy, ref_destruct = get_inputs(ref_path)
    lumpy, destruct = get_inputs(ref_path)

    assert os.path.exists(ref_lumpy)
    assert os.path.exists(lumpy)
    # assert os.path.exists(ref_destruct)
    # assert os.path.exists(destruct)

    compare_breakpoint_calling(ref_lumpy, lumpy, ref_destruct, destruct)
