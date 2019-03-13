'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import csv
import gzip
import yaml
import logging
import pandas as pd
from single_cell.utils import helpers
import time


def generate_dtype_yaml(csv_file, yaml_filename):
    pandas_to_std_types = {
        "bool": "bool",
        "int64": "int",
        "float64": "float",
        "object": "str",
    }

    if isinstance(csv_file, str):
        print csv_file
        print helpers.get_compression_type_pandas(csv_file)
        print os.stat(csv_file)
        data = pd.read_csv(
            csv_file, compression=helpers.get_compression_type_pandas(csv_file)
        )
    elif isinstance(csv_file, pd.DataFrame):
        data = csv_file
    else:
        raise ValueError(
            "Incorrect input data type. must be a dataframe or csv file path"
        )

    if len(data.columns) != len(data.columns.unique()):
        raise ValueError('duplicate columns not supported')

    typeinfo = {}
    for column, dtype in data.dtypes.iteritems():
        typeinfo[column] = pandas_to_std_types[str(dtype)]

    with open(yaml_filename, 'w') as f:
        yaml.dump(typeinfo, f, default_flow_style=False)


def annotate_metrics(infile, sample_info, outfile, yamlfile=None):
    metrics_df = pd.read_csv(infile)

    cells = metrics_df["cell_id"]

    for cell in cells:
        coldata = sample_info[cell]

        for column, value in coldata.iteritems():
            metrics_df.loc[metrics_df["cell_id"] == cell, column] = value

    write_to_file(
        metrics_df, outfile, index=False,
        compression=helpers.get_compression_type_pandas(outfile)
    )

    if yamlfile:
        generate_dtype_yaml(metrics_df, yamlfile)


def concatenate_csv(in_filenames, out_filename, nan_val='NA', key_column=None, sep=',', yamlfile=None):
    data = []

    if not isinstance(in_filenames, dict):
        in_filenames = dict(enumerate(in_filenames))

    for key, in_filename in in_filenames.iteritems():
        with open(in_filename) as f:
            first_line = f.readline()
            if len(first_line) == 0:
                continue
        df = pd.read_csv(in_filename, dtype=str, sep=sep)
        if key_column is not None:
            df[key_column] = str(key)
        data.append(df)
    data = pd.concat(data, ignore_index=True)
    data = data.fillna(nan_val)

    write_to_file(data, out_filename, index=False, compression=helpers.get_compression_type_pandas(out_filename))

    if yamlfile:
        generate_dtype_yaml(data, yamlfile)



def merge_csv(in_filenames, out_filename, how, on, nan_val='NA', sep=',', suffixes=None, yamlfile=None):
    data = []

    if isinstance(in_filenames, dict):
        in_filenames = in_filenames.values()

    for in_filename in in_filenames:
        print 'merging', in_filename
        with open(in_filename) as f:
            first_line = f.readline()
            if len(first_line) == 0:
                continue
        data.append(pd.read_csv(in_filename, sep=sep, dtype=str))

    data = merge_frames(data, how, on, suffixes = suffixes)
    data = data.fillna(nan_val)

    write_to_file(data, out_filename, index=False, compression=helpers.get_compression_type_pandas(out_filename))

    if yamlfile:
        generate_dtype_yaml(data, yamlfile)



def write_to_file(df, out_filename, index=True, compression=None, sep=','):
    if compression == "hdf5":
        df.to_hdf5(out_filename, index=index)
    else:
        df.to_csv(out_filename, na_rep='NA', index=index, compression=helpers.get_compression_type_pandas(out_filename), sep=sep)


def merge_frames(frames, how, on, suffixes=None):
    '''
    annotates input_df using ref_df
    '''

    suff = ['','']

    if ',' in on:
        on = on.split(',')

    if len(frames) == 1:
        return frames[0]
    else:
        left = frames[0]
        right = frames[1]

        if suffixes:
            suff = (suffixes[0],suffixes[1])

        merged_frame = pd.merge(left, right,
                                how=how,
                                on=on,
                                suffixes=suff)
        for i,frame in enumerate(frames[2:]):

            if suffixes:
                suff = (suffixes[i+2],suffixes[i+2])

            merged_frame = pd.merge(merged_frame, frame,
                                    how=how,
                                    on=on,
                                    suffixes=suff)
        return merged_frame


def concatenate_csv_lowmem(in_filenames, out_filename, yamlfile=None):
    """merge csv files, uses csv module to handle inconsistencies in column
    indexes, pandas uses a lot of memory
    :param in_filenames: input file dict
    :param out_filename: output file
    """
    if isinstance(in_filenames, dict):
        in_filenames = in_filenames.values()

    compression = helpers.get_file_format(out_filename)
    if compression == "gzip":
        out_writer = gzip.open(out_filename, 'w')
    else:
        out_writer = open(out_filename, 'w')

    writer = None


    for infile in in_filenames:
        if helpers.is_gzip(infile):
            inp = gzip.open(infile)
        else:
            inp = open(infile)

        reader = csv.DictReader(inp)

        for row in reader:
            if not writer:
                if compression == "gzip":
                    writer = csv.DictWriter(out_writer,
                                            fieldnames=reader._fieldnames)
                else:
                    writer = csv.DictWriter(out_writer,
                                            fieldnames=reader._fieldnames)
                writer.writeheader()

            print row
            writer.writerow(row)
        inp.close()

    if not writer:
        logging.getLogger("single_cell.helpers.csvutils").warn(
            "no data to merge, generating an empty file"
        )

        #if inputs have headers write header to output
        if reader._fieldnames:
            writer = csv.DictWriter(out_writer,
                                    fieldnames=reader._fieldnames)
            writer.writeheader()

    out_writer.close()

    if yamlfile:
        generate_dtype_yaml(out_filename, yamlfile)

