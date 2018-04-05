'''
Created on Jul 24, 2017

@author: dgrewal
'''
import pandas as pd
from scripts import MergeFiles
from scripts import GetCounts

from single_cell.utils import csvutils

def concat_csv(in_filenames, out_filename):
    csvutils.concatenate_csv_lowmem(in_filenames, out_filename)


def merge_csv(in_filenames, out_filename, how, on, nan_val='NA', sep=',', suffixes=None):
    csvutils.merge_csv(in_filenames, out_filename, how, on, sep=sep, suffixes=suffixes)


def get_counts(bam, bai, positions, output, cell_id):
    counts = GetCounts(bam, positions, output, cell_id)
    counts.main()
