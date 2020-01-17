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
    
    museq = os.path.join(path, lib + "snv_museq.csv.gz") 

    strelka = os.path.join(path, lib + "snv_strelka.csv.gz")

    snpeff = os.path.join(path, lib + "snv_snpeff.csv.gz")

    return museq, strelka, snpeff

def compare_variant_calling(ref_museq, museq, 
                            ref_strelka, strelka,
                            ref_snpeff, snpeff):

    compare.compare_variant_calls(ref_museq, museq)
    compare.compare_variant_calls(ref_strelka, strelka)
    compare.compare_annotation(ref_snpeff, snpeff)

if __name__ == "__main__":

    output_path = sys.argv[1]
    ref_path = sys.argv[2]

    ref_museq, ref_strelka, ref_snpeff = get_inputs(ref_path)
    museq, strelka, snpeff = get_inputs(output_path)

    compare_variant_calling(ref_museq, museq,
                    ref_strelka, strelka,
                    ref_snpeff, snpeff)
