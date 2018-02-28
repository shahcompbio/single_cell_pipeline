'''
Created on Jul 24, 2017

@author: dgrewal
'''
import pandas as pd
from scripts import MergeFiles
from scripts import GetCounts

from single_cell.utils import csvutils

def merge_csv(in_filenames, out_filename, how, on, nan_val='NA', sep=','):
    csvutils.merge_csv(in_filenames, out_filename, how, on, sep=sep)

def get_counts(bam, positions, output, sample_id):
    counts = GetCounts(bam, positions, output, sample_id)
    counts.main()
