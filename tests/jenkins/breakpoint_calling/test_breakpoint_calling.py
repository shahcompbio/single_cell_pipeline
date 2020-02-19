import os.path
import sys
from single_cell.utils import compare
from single_cell.utils import csvutils

def get_inputs(path):
    """"
    get metrics and gc metrics given a directory and library
    :param path:  path to metrics files
    """

    must_exist = ["destruct_breakpoints_library.csv.gz",
                  "destruct_breakpoints_library.csv.gz.yaml",
                  "destruct_cell_counts.csv.gz",
                  "destruct_cell_counts.csv.gz.yaml",
                  "input.yaml",
                  "lumpy_breakpoints.bed",
                  "lumpy_breakpoints_evidence.csv.gz",
                  "lumpy_breakpoints_evidence.csv.gz.yaml",
                  "metadata.yaml"]

    lumpy_breakpoints = os.path.join(path, "lumpy_breakpoints.csv.gz")
    destruct_breakpoints = os.path.join(path, "destruct_breakpoints.csv.gz")

    must_exist = [os.path.join(path, f) for f in must_exist]

    return must_exist, lumpy_breakpoints, destruct_breakpoints


def test_breakpoint_calling(args):
    output_path = args[1]
    ref_path = args[2]

    ref_must_exist, ref_lumpy, ref_destruct = get_inputs(ref_path)
    must_exist, lumpy, destruct = get_inputs(output_path)

    assert all(map(os.path.exists, ref_must_exist))
    assert all(map(os.path.exists, must_exist))

    compare.compare_breakpoint_calls(ref_lumpy, lumpy)

    ref_destruct = csvutils.read_csv_and_yaml(ref_destruct)
    destruct = csvutils.read_csv_and_yaml(destruct)

    assert ref_destruct.empty and destruct.empty

if __name__ == "__main__":
    test_breakpoint_calling(sys.argv)
