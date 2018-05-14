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

    biowrappers_hdf5.concatenate_tables(
        in_files,
        out_file,
        drop_duplicates=drop_duplicates,
        in_memory=in_memory,
        non_numeric_as_category=non_numeric_as_category)


def merge_cells_in_memory(hdf_input, output_store_obj, tablename, dtypes={}):
    data = []

    with pd.HDFStore(hdf_input, 'r') as input_store:
        for tableid in input_store:
            data.append(input_store[tableid])

    data = pd.concat(data)
    data = data.reset_index()

    for col, dtype in dtypes.iteritems():
        data[col] = data[col].astype(dtype)

    output_store_obj.put(tablename, data, format="table")


def merge_cells_on_disk(hdf_input, output_store_obj, tablename, dtypes={}):

    with pd.HDFStore(hdf_input, 'r') as input_store:
        for tableid in input_store:
            celldata = input_store[tableid]

            for col, dtype in dtypes.iteritems():
                celldata[col] = celldata[col].astype(dtype)

            if tablename not in output_store_obj:
                output_store_obj.put(tablename, celldata, format='table')
            else:
                output_store_obj.append(tablename, celldata, format='table')


def merge_per_cell_tables(
        infile, output, tablename, in_memory=True, dtypes={}):

    if isinstance(output, pd.HDFStore):
        output_store = output
    else:
        output_store = pd.HDFStore(output, 'w', complevel=9, complib='blosc')

    if in_memory:
        merge_cells_in_memory(infile, output_store, tablename, dtypes=dtypes)
    else:
        merge_cells_on_disk(
            infile,
            output_store,
            tablename,
            dtypes=dtypes)

    if not isinstance(output, pd.HDFStore):
        output_store.close()


def annotate_per_cell_store_with_dict(infile, annotation_data, output):
    """
    adds new cols to dataframes in store from a dictionary
    store must be split by cells (no tables with all cells merged together)
    annotation_data must be 2 level dict with cellid as first key,
    colnames as second, values for cols as leaves

    """
    with pd.HDFStore(output, 'w', complevel=9, complib='blosc') as output, pd.HDFStore(infile) as input_store:
        for tableid in input_store.keys():
            data = input_store[tableid]

            cell_id = data["cell_id"].iloc[0]

            cell_info = annotation_data[cell_id]

            for colname, value in cell_info.iteritems():
                data[colname] = value

                output.put(tableid, data, format="table")


def annotate_store_with_dict(infile, annotation_data, output):
    """
    adds new cols to dataframes in store from a dictionary
    """
    with pd.HDFStore(output, 'w', complevel=9, complib='blosc') as output, pd.HDFStore(infile) as input_store:
        for tableid in input_store.keys():
            data = input_store[tableid]

            cells = data["cell_id"].unique()

            for cellid in cells:
                cell_info = annotation_data[cellid]
                for colname, value in cell_info.iteritems():
                    data.loc[data["cell_id"] == cellid, colname] = value

            output.put(tableid, data, format="table")
