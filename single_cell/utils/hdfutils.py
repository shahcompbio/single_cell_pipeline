'''
Created on May 9, 2018

@author: dgrewal
'''
import gc
import logging

import pandas as pd
from biowrappers.components.io.hdf5 import tasks as biowrappers_hdf5
from single_cell.utils import csvutils
from single_cell.utils import helpers


def get_min_itemsize(files):
    min_itemsize = {}

    itemsizes = biowrappers_hdf5._get_min_itemsize(files)

    for tablename, tableitemsizes in itemsizes.items():
        for col, len_item in tableitemsizes.items():
            if col in min_itemsize:
                min_itemsize[col] = max(min_itemsize[col], len_item)
            else:
                min_itemsize[col] = len_item

            min_itemsize[col] += 2
    return min_itemsize


def cast_columns(df):
    if not isinstance(df, pd.DataFrame):
        return df

    reference, ignore_cols = helpers.get_coltype_reference()

    for col in df.columns.values:
        if col in reference:
            try:
                df[col] = df[col].astype(reference[col])
            except ValueError as exc:
                logging.getLogger("single_cell.helpers.hdfutils").warn(
                    "could not cast {} due to error: {}".format(col, exc.message)
                )
        else:
            if col in ignore_cols:
                logging.getLogger("single_cell.helpers.hdfutils").warn(
                    'Could not cast {}, please add the expected data type to reference'.format(col)
                )

    return df


def cast_h5_file(h5data, output):
    with pd.HDFStore(output, 'w', complevel=9, complib='blosc') as output:
        with pd.HDFStore(h5data, 'r') as input:
            tablenames = input.keys()
            for tablename in tablenames:
                df = input[tablename]
                df = cast_columns(df)
                output.put(tablename, df, format='table')
                # helps lower mem usage
                df = None
                gc.collect()


def concat_csvs_to_hdf(infiles, outfile, tablenames):
    with pd.HDFStore(outfile, 'w', complevel=9, complib='blosc') as output:
        for infile, tablename in zip(infiles, tablenames):
            df = pd.read_csv(infile)
            # pytables silently ignores empty tables
            if df.empty:
                df = df.append(pd.Series([]), ignore_index=True)
            output.put(tablename, df, format='table')


def merge_csvs_to_hdf_in_memory(infiles, outfile, tablename):
    data = [pd.read_csv(infile) for infile in infiles]
    data = pd.concat(data)

    with pd.HDFStore(outfile, 'w', complevel=9, complib='blosc') as output:
        output.put(tablename, data, format='table')


def merge_csvs_to_hdf_on_disk(infiles, outfile, tablename):
    with pd.HDFStore(outfile, 'w', complevel=9, complib='blosc') as output:

        for infile in infiles:
            data = pd.read_csv(infile)

            if tablename not in output:
                output.put(tablename, data, format='table')
            else:
                output.append(tablename, data, format='table')


def convert_csv_to_hdf(infile, outfile, tablename):
    df = pd.read_csv(infile)
    df = df.infer_objects()

    with pd.HDFStore(outfile, 'w', complevel=9, complib='blosc') as out_store:
        out_store.put(tablename, df, format='table')


def set_categories_df(df, categories):
    if not isinstance(df, pd.DataFrame):
        return df
    for colname in df.columns.values:
        if colname in categories:
            df[colname] = df[colname].cat.set_categories(categories[colname])
    return df


def concat_hdf_tables(in_files, out_file, categories={}):
    chunksize = 10 ** 6

    with pd.HDFStore(out_file, 'w', complevel=9, complib='blosc') as output:
        for infile in in_files:
            with pd.HDFStore(infile, 'r') as input_store:
                tables = input_store.keys()

            for table in tables:
                # metadata columns from original
                # files can cause loading errors
                if table.endswith("meta"):
                    continue

                for chunk in pd.read_hdf(
                        infile, key=table, chunksize=chunksize):

                    chunk = set_categories_df(chunk, categories)

                    if table not in output:
                        output.put(table, chunk, format='table')
                    else:
                        output.append(table, chunk, format='table')


def merge_cells_in_memory(
        hdf_input, output_store_obj, tablename,
        tables_to_merge=None, dtypes={}):
    data = []

    with pd.HDFStore(hdf_input, 'r') as input_store:
        if not tables_to_merge:
            tables_to_merge = input_store.keys()

        for tableid in tables_to_merge:
            data.append(input_store[tableid])

    data = pd.concat(data)
    data = data.reset_index()

    for col, dtype in dtypes.items():
        data[col] = data[col].astype(dtype)

    output_store_obj.put(tablename, data, format="table")


def merge_cells_on_disk(
        hdf_input, output_store_obj, tablename,
        tables_to_merge=None, dtypes={}):
    with pd.HDFStore(hdf_input, 'r') as input_store:

        if not tables_to_merge:
            tables_to_merge = input_store.keys()

        for tableid in tables_to_merge:
            celldata = input_store[tableid]

            for col, dtype in dtypes.items():
                celldata[col] = celldata[col].astype(dtype)

            if tablename not in output_store_obj:
                output_store_obj.put(tablename, celldata, format='table')
            else:
                output_store_obj.append(tablename, celldata, format='table')


def merge_per_cell_tables(
        infile, output, out_tablename,
        tables_to_merge=None, in_memory=True, dtypes={}):
    if isinstance(output, pd.HDFStore):
        output_store = output
    else:
        output_store = pd.HDFStore(output, 'w', complevel=9, complib='blosc')

    if in_memory:
        merge_cells_in_memory(
            infile,
            output_store,
            out_tablename,
            tables_to_merge=tables_to_merge,
            dtypes=dtypes)
    else:
        merge_cells_on_disk(
            infile,
            output_store,
            out_tablename,
            tables_to_merge=tables_to_merge,
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

            for colname, value in cell_info.items():
                data[colname] = value

                output.put(tableid, data, format="table")


def annotate_store_with_dict(infile, annotation_data, output, tables=None):
    """
    adds new cols to dataframes in store from a dictionary
    """

    if isinstance(output, pd.HDFStore):
        output_store = output
    else:
        output_store = pd.HDFStore(output, 'w', complevel=9, complib='blosc')

    if isinstance(infile, pd.HDFStore):
        input_store = infile
    else:
        input_store = pd.HDFStore(infile, 'r')

    if not tables:
        tables = input_store.keys()

    for tableid in tables:
        data = input_store[tableid]

        cells = data["cell_id"].unique()

        for cellid in cells:
            cell_info = annotation_data[cellid]
            for colname, value in cell_info.items():
                data.loc[data["cell_id"] == cellid, colname] = value

        output_store.put(tableid, data, format="table")

    if not isinstance(output, pd.HDFStore):
        output_store.close()

    if not isinstance(input, pd.HDFStore):
        input_store.close()


def convert_hdf_to_csv(h5_input, outputs, dtypes, chunksize=10**6):
    for tablename, outfile in outputs.items():
        header=False
        for chunk in pd.read_hdf(h5_input, key=tablename, chunksize=chunksize):
            compression = helpers.get_compression_type_pandas(outfile)

            if header:
                chunk.to_csv(outfile, index=False, header=False, mode='a', compression=compression)
            else:
                chunk.to_csv(outfile, index=False, mode='w', compression=compression)
                header = True

        csvutils.write_metadata(outfile, dtypes)
