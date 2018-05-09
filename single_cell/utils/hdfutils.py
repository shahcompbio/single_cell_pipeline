'''
Created on May 9, 2018

@author: dgrewal
'''
import pandas as pd

from biowrappers.components.io.hdf5 import tasks as biowrappers_hdf5



def convert_csv_to_hdf(infile, outfile, tablename):
    df = pd.read_csv(infile)

    df = df.infer_objects()

    with pd.HDFStore(outfile, 'w', complevel=9, complib='blosc') as out_store:
        out_store.put(tablename, df, format='table')


def concat_hdf_tables(in_files, out_file, drop_duplicates=False,
                     in_memory=True, non_numeric_as_category=True):

    biowrappers_hdf5.concatenate_tables(in_files, out_file, drop_duplicates=drop_duplicates, in_memory=in_memory, non_numeric_as_category=non_numeric_as_category)
