'''
Created on Feb 19, 2018

@author: dgrewal
'''

import pandas as pd
import csv
import warnings


def concatenate_csv(in_filenames, out_filename, nan_val='NA'):
    data = []

    if isinstance(in_filenames, dict):
        in_filenames = in_filenames.values()

    
    for in_filename in in_filenames:
        with open(in_filename) as f:
            first_line = f.readline()
            if len(first_line) == 0:
                continue
        data.append(pd.read_csv(in_filename, dtype=str))
    data = pd.concat(data, ignore_index=True)
    data = data.fillna(nan_val)
    data.to_csv(out_filename, index=False)


def merge_csv(in_filenames, out_filename, how, on, nan_val='NA', sep=','):
    data = []

    if isinstance(in_filenames, dict):
        in_filenames = in_filenames.values()

    for in_filename in in_filenames:
        with open(in_filename) as f:
            first_line = f.readline()
            if len(first_line) == 0:
                continue
        data.append(pd.read_csv(in_filename, sep=sep, dtype=str))

    data = merge_frames(data, how, on)
    data = data.fillna(nan_val)
    data.to_csv(out_filename, index=False)


def merge_frames(frames, how, on):
    '''
    annotates input_df using ref_df
    '''

    if ',' in on:
        on = on.split(',')

    if len(frames) == 1:
        return frames[0]
    else:
        left = frames[0]
        right = frames[1]
        merged_frame = pd.merge(left, right,
                                how=how,
                                on=on)
        for frame in frames[2:]:
            merged_frame = pd.merge(merged_frame, frame,
                                    how=how,
                                    on=on)
        return merged_frame


def concatenate_csv_lowmem(in_filenames, out_filename):
    """merge csv files, uses csv module to handle inconsistencies in column
    indexes, pandas uses a lot of memory
    :param in_filenames: input file dict
    :param out_filename: output file
    """
    writer = None
    for _,infile in in_filenames.iteritems():

        with open(infile) as inp:
            reader= csv.DictReader(inp)

            for row in reader:
                if not writer:
                    writer = csv.DictWriter(open(out_filename, "w"),
                                            fieldnames=reader._fieldnames)
                    writer.writeheader()

                writer.writerow(row)

    if not writer:
        warnings.warn("no data to merge, generating an empty file")
        open(out_filename, 'w').close()
