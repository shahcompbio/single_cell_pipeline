'''
Created on Feb 19, 2018

@author: dgrewal
'''

import pandas as pd
import csv
import warnings
import gzip

def annotate_metrics(infile, sample_info, outfile, compression=None):
    metrics_df = pd.read_csv(infile)

    cols = ["pick_met", "condition", "sample_type",
            "index_i5", "index_i7", "row", "column",
            'img_col', "primer_i5", "primer_i7",
            ]

    for col in cols:
        metrics_df[col] = "NA"

    cells = metrics_df["cell_id"]

    for col in cols:
        coldata = [sample_info[cell][col] for cell in cells]
        metrics_df[col] = coldata


    write_to_file(metrics_df, outfile, index=False, compression=compression)


def concatenate_csv(in_filenames, out_filename, nan_val='NA', compression=None, key_column=None, sep=','):
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

    write_to_file(data, out_filename, index=False, compression=compression)


def merge_csv(in_filenames, out_filename, how, on, nan_val='NA', sep=',', suffixes=None, compression=None):
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

    write_to_file(data, out_filename, index=False, compression=compression)


def write_to_file(df, out_filename, index=True, compression=None, sep=','):
    if compression == "hdf5":
        df.to_hdf5(out_filename, index=index)
    else:
        df.to_csv(out_filename, index=index, compression=compression, sep=sep)


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


def concatenate_csv_lowmem(in_filenames, out_filename, compression=None):
    """merge csv files, uses csv module to handle inconsistencies in column
    indexes, pandas uses a lot of memory
    :param in_filenames: input file dict
    :param out_filename: output file
    """
    if compression and not compression == "gzip":
        raise Exception("doesn't support {} compression. only supports gzip".format(compression))

    writer = None

    if isinstance(in_filenames, dict):
        in_filenames = in_filenames.values()

    for infile in in_filenames:

        with open(infile) as inp:
            reader= csv.DictReader(inp)

            for row in reader:
                if not writer:
                    if compression == "gzip":
                        writer = csv.DictWriter(gzip.open(out_filename, "w"),
                                                fieldnames=reader._fieldnames)
                    else:
                        writer = csv.DictWriter(open(out_filename, "w"),
                                                fieldnames=reader._fieldnames)
                    writer.writeheader()

                writer.writerow(row)

    if not writer:
        warnings.warn("no data to merge, generating an empty file")

        #if inputs have headers write header to output
        if reader._fieldnames:
            writer = csv.DictWriter(open(out_filename, "w"),
                                    fieldnames=reader._fieldnames)
            writer.writeheader()
        else:
            open(out_filename, 'w').close()
