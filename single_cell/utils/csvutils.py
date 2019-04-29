'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import shutil

import pandas as pd
import yaml
from single_cell.utils import helpers


def read_csv_and_yaml_by_chunks(infile, dtypes, columns, chunksize, header):
    dfs = pd.read_csv(
        infile, compression=helpers.get_compression_type_pandas(infile),
        dtype=dtypes, header=header, chunksize=chunksize
    )
    for data in dfs:
        if header is None:
            data.columns = columns
        else:
            assert list(data.columns.values) == columns
        yield data


def read_csv_and_yaml(infile, chunksize=None):
    with open(infile) as f:
        first_line = f.readline()
        if len(first_line) == 0:
            return

    header, dtypes, columns = get_metadata(infile)

    # if header exists then use first line (0) as header
    header = 0 if header else None

    if chunksize:
        return read_csv_and_yaml_by_chunks(infile, dtypes, columns, chunksize, header)
    else:
        try:
            data = pd.read_csv(
                infile, compression=helpers.get_compression_type_pandas(infile),
                dtype=dtypes, header=header
            )
        except pd.errors.EmptyDataError:
            data = pd.DataFrame(columns=columns)
        if header is None:
            data.columns = columns
        else:
            assert list(data.columns.values) == columns
        return data


def get_metadata(filepath):
    if not os.path.exists(filepath + '.yaml'):
        with helpers.getFileHandle(filepath) as inputfile:
            columns = inputfile.readline().strip().split(',')
            header = True
            dtypes = None
            return header, dtypes, columns

    with open(filepath + '.yaml') as yamlfile:
        yamldata = yaml.load(yamlfile)

    header = yamldata['header']

    dtypes = {}
    columns = []
    for coldata in yamldata['columns']:
        colname = coldata['name']

        dtypes[colname] = coldata['dtype']

        columns.append(colname)

    return header, dtypes, columns


def load_csv_metadata(csvfile):
    yamlfile = csvfile + '.yaml'

    if not os.path.exists(yamlfile):
        return None

    with open(yamlfile) as yamlinput:
        return yaml.load(yamlinput)


def write_dataframe_to_csv_and_yaml(df, outfile):
    compression = helpers.get_compression_type_pandas(outfile)

    if compression == 'h5':
        df.to_hdf5(outfile)
    else:
        df.to_csv(outfile, compression=compression, header=False, na_rep='NA', index=False)
        generate_yaml_for_csv(df, outfile + '.yaml')


def generate_yaml_for_csv(filepath, outputyaml, header=False, columns=None):
    if isinstance(filepath, pd.DataFrame):
        types = generate_dtype_yaml(filepath)
        columns = list(filepath.columns.values)
    else:
        with helpers.getFileHandle(filepath) as infile:
            if not columns:
                columns = infile.readline().strip().split(',')
            types = generate_dtype_yaml(filepath, columns=columns)

    write_to_yaml(header, types, columns, outputyaml)


def write_to_yaml(header, types, columns, outputyaml):
    yamldata = {'header': header, 'columns': []}

    for column in columns:
        data = {'name': column, 'dtype': types[column]}
        yamldata['columns'].append(data)

    with open(outputyaml, 'w') as f:
        yaml.dump(yamldata, f, default_flow_style=False)


def generate_dtype_yaml(csv_file, yaml_filename=None, columns=None):
    pandas_to_std_types = {
        "bool": "bool",
        "int64": "int",
        "float64": "float",
        "object": "str",
    }

    if isinstance(csv_file, str):
        # dont need to read entire file to generate
        # this yaml. read the first chunk only
        chunksize = 10 ** 6
        data = pd.read_csv(
            csv_file, compression=helpers.get_compression_type_pandas(csv_file),
            chunksize=chunksize
        )
        data = next(data)
        if columns:
            data.columns = columns
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
        if column in ['chr', 'chrom', 'chromosome']:
            typeinfo[column] = 'str'
        else:
            typeinfo[column] = pandas_to_std_types[str(dtype)]

    if yaml_filename:
        with open(yaml_filename, 'w') as f:
            yaml.dump(typeinfo, f, default_flow_style=False)
    else:
        return typeinfo


def annotate_metrics(infile, sample_info, outfile, yamlfile=None):
    metrics_df = read_csv_and_yaml(infile)

    cells = metrics_df["cell_id"]

    for cell in cells:
        coldata = sample_info[cell]

        for column, value in coldata.iteritems():
            metrics_df.loc[metrics_df["cell_id"] == cell, column] = value

    write_dataframe_to_csv_and_yaml(metrics_df, outfile)


def concatenate_csv(in_filenames, out_filename, nan_val='NA', key_column=None):
    data = []

    if not isinstance(in_filenames, dict):
        in_filenames = dict(enumerate(in_filenames))

    for key, in_filename in in_filenames.iteritems():
        with open(in_filename) as f:
            first_line = f.readline()
            if len(first_line) == 0:
                continue
        df = read_csv_and_yaml(in_filename)
        if key_column is not None:
            df[key_column] = str(key)
        data.append(df)
    data = pd.concat(data, ignore_index=True)
    data = data.fillna(nan_val)

    write_dataframe_to_csv_and_yaml(data, out_filename)


def extrapolate_types_from_yaml_files(csv_files, outputyaml):
    precedence = ['str', 'float', 'int', 'bool']

    yamldata = [get_metadata(csvfile) for csvfile in csv_files]

    header = set([val[0] for val in yamldata])
    assert len(header) == 1, 'mimatched yaml files'
    header = list(header)[0]

    cols = set([tuple(val[2]) for val in yamldata])
    assert len(cols) == 1, 'mimatched yaml files'
    cols = list(cols)[0]

    dtypes = [val[1] for val in yamldata]

    merged_dtypes = {}
    for dtype_val in dtypes:
        for col, dtype in dtype_val.items():
            if col not in merged_dtypes:
                merged_dtypes[col] = dtype
            else:
                og_index = precedence.index(merged_dtypes[col])
                new_index = precedence.index(dtype)
                if new_index < og_index:
                    merged_dtypes[col] = dtype

    write_to_yaml(header, merged_dtypes, cols, outputyaml)


def concatenate_csv_files_quick_lowmem(inputfiles, output):
    if isinstance(inputfiles, dict):
        inputfiles = inputfiles.values()

    with helpers.getFileHandle(output, 'w') as outfile:
        for infile in inputfiles:
            with helpers.getFileHandle(infile) as inputdata:
                shutil.copyfileobj(inputdata, outfile, length=16 * 1024 * 1024)

    extrapolate_types_from_yaml_files(inputfiles, output+'.yaml')


def merge_csv(in_filenames, out_filename, how, on, nan_val='NA', suffixes=None):
    data = []

    if isinstance(in_filenames, dict):
        in_filenames = in_filenames.values()

    for in_filename in in_filenames:
        indata = read_csv_and_yaml(in_filename)
        if indata:
            data.append(indata)

    data = merge_frames(data, how, on, suffixes=suffixes)
    data = data.fillna(nan_val)

    write_dataframe_to_csv_and_yaml(data, out_filename)


def merge_frames(frames, how, on, suffixes=None):
    '''
    annotates input_df using ref_df
    '''

    suff = ['', '']

    if ',' in on:
        on = on.split(',')

    if len(frames) == 1:
        return frames[0]
    else:
        left = frames[0]
        right = frames[1]

        if suffixes:
            suff = (suffixes[0], suffixes[1])

        merged_frame = pd.merge(left, right,
                                how=how,
                                on=on,
                                suffixes=suff)
        for i, frame in enumerate(frames[2:]):

            if suffixes:
                suff = (suffixes[i + 2], suffixes[i + 2])

            merged_frame = pd.merge(merged_frame, frame,
                                    how=how,
                                    on=on,
                                    suffixes=suff)
        return merged_frame


def finalize_csv(infile, outfile):
    header, dtypes, columns = get_metadata(infile)

    assert header==False, 'file already contains a header'

    header = ','.join(columns) + '\n'

    with helpers.getFileHandle(outfile, 'w') as output:
        output.write(header)
        with helpers.getFileHandle(infile) as indata:
            shutil.copyfileobj(indata, output, length=16 * 0124 * 1024)

    write_to_yaml(True, dtypes, columns, outfile+'.yaml')


def prep_csv_files(filepath, outputfile):
    outputyaml = outputfile + '.yaml'
    with helpers.getFileHandle(outputfile, 'w') as out_writer:
        generate_yaml_for_csv(filepath, outputyaml)

        with helpers.getFileHandle(filepath) as infile:
            # skip header
            infile.readline()
            shutil.copyfileobj(infile, out_writer, length=16 * 1024 * 1024)
