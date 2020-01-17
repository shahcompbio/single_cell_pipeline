import sys 
from single_cell.utils import compare

def compare_merge_wgs(output_path, ref_counts):
    compare.assess_bam_files(output_path, ref_counts)
    
if __name__ == "__main__":

    output_path = sys.argv[1]
    ref_counts = sys.argv[2]

    compare_merge_wgs(output_path, ref_counts)
