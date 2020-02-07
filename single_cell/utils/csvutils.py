import gzip
import os
import shutil
import itertools
import pandas as pd
import yaml


class CsvMergeDtypesEmptyMergeSet(Exception):
    pass

class CsvConcatNaNIntDtypeException(Exception):
    pass


class CsvMergeCommonColException(Exception):
    pass


class DtypesMergeException(Exception):
    pass


class CsvConcatException(Exception):
    pass


class CsvAnnotateError(Exception):
    pass


class CsvMergeException(Exception):
    pass


class CsvMergeColumnMismatchException(Exception):
    pass


class CsvParseError(Exception):
    pass


class CsvInputError(Exception):
    pass


def pandas_to_std_types():
    return {
        "bool": "bool",
        "int64": "int",
        "int": "int",
        "Int64": "int",
        "float64": "float",
        "float": "float",
        "object": "str",
        "str": "str",
        "category": "str",
        "NaN": "NaN",
    }


class CsvWriterError(Exception):
    pass


class CsvTypeMismatch(Exception):
    def __init__(self, column, expected_dtype, dtype):
        self.column = column
        self.dtype = dtype
        self.expected_dtype = expected_dtype

    def __str__(self):
        message = 'mismatching types for col {}. types were {} and {}'
        message = message.format(self.column, self.expected_dtype, self.dtype)
        return message


class IrregularCsvInput(object):
    def __init__(self, filepath, na_rep='NaN'):
        """
        csv file and all related metadata
        :param filepath: path to csv
        :type filepath: str
        :param na_rep: replace na with this
        :type na_rep: str
        """
        self.filepath = filepath

        self.na_rep = na_rep

        metadata = self.__generate_metadata()

        self.header, self.sep, self.dtypes, self.columns = metadata

        self.__confirm_compression_type_pandas()

    @property
    def yaml_file(self):
        return self.filepath + '.yaml'

    def get_dtypes_from_df(self, df):
        type_converter = pandas_to_std_types()

        typeinfo = {}
        for column, dtype in df.dtypes.items():
            if column in ['chr', 'chrom', 'chromosome']:
                typeinfo[column] = 'str'
                typeinfo[column] = 'str'
            else:
                if df.empty:
                    typeinfo[column] = self.na_rep
                else:
                    typeinfo[column] = type_converter[str(dtype)]

        return typeinfo

    def __detect_sep_from_header(self, header):
        """
        detect whether file is tab or comma separated from header
        :param header: header line
        :type header: str
        :return: separator
        :rtype: str
        """
        if '\t' in header and ',' in header:
            raise CsvParseError("Unable to detect separator from {}".format(header))

        if '\t' not in header and ',' not in header:
            raise CsvParseError("Unable to detect separator from {}".format(header))

        if '\t' in header:
            return '\t'
        elif ',' in header:
            return ','

    def __confirm_compression_type_pandas(self):
        filepath = self.filepath
        if filepath.endswith('.tmp'):
            filepath = filepath[:-4]

        _, ext = os.path.splitext(filepath)

        if not ext == ".gz":
            raise CsvInputError("{} is not supported".format(ext))

    def __generate_metadata(self):
        with gzip.open(self.filepath, 'rt') as inputfile:
            header = inputfile.readline().strip()
            sep = self.__detect_sep_from_header(header)
            columns = header.split(sep)
            header = True
            dtypes = self.__generate_dtypes(sep=sep)
            return header, sep, dtypes, columns

    def __generate_dtypes(self, columns=None, sep=','):
        data = pd.read_csv(
            self.filepath, compression='gzip', chunksize=10 ** 6,
            sep=sep
        )
        data = next(data)

        if columns:
            data.columns = columns

        typeinfo = self.get_dtypes_from_df(data)
        return typeinfo


class CsvInput(object):
    def __init__(self, filepath, na_rep='NaN'):
        """
        csv file and all related metadata
        :param filepath: path to csv
        :type filepath: str
        :param na_rep: replace na with this
        :type na_rep: str
        """
        self.filepath = filepath

        self.na_rep = na_rep

        self.sep = ','

        metadata = self.__parse_metadata()

        self.header, self.dtypes, self.columns = metadata

        self.__confirm_compression_type_pandas()

    def cast_dataframe(self, df):
        for column_name in df.columns.values:
            dtype = self.dtypes[column_name]
            df[column_name] = df[column_name].astype(dtype)
        return df

    @property
    def yaml_file(self):
        return self.filepath + '.yaml'

    def __confirm_compression_type_pandas(self):
        filepath = self.filepath
        if filepath.endswith('.tmp'):
            filepath = filepath[:-4]

        _, ext = os.path.splitext(filepath)

        if not ext == ".gz":
            raise CsvInputError("{} is not supported".format(ext))

    def __parse_metadata(self):
        with open(self.filepath + '.yaml') as yamlfile:
            yamldata = yaml.safe_load(yamlfile)

        header = yamldata['header']
        sep = yamldata['sep']

        assert sep == self.sep

        dtypes = {}
        columns = []
        for coldata in yamldata['columns']:
            colname = coldata['name']

            dtypes[colname] = coldata['dtype']

            columns.append(colname)

        return header, dtypes, columns

    def __verify_data(self, df):
        if not set(list(df.columns.values)) == set(self.columns):
            raise CsvParseError("metadata mismatch in {}".format(self.filepath))

    def read_csv(self, chunksize=None):
        def return_gen(df_iterator):
            for df in df_iterator:
                self.__verify_data(df)
                yield df

        dtypes = {k: v for k, v in self.dtypes.items() if v != "NA"}
        # if header exists then use first line (0) as header
        header = 0 if self.header else None
        names = None if self.header else self.columns

        try:
            data = pd.read_csv(
                self.filepath, compression='gzip', chunksize=chunksize,
                sep=self.sep, header=header, names=names, dtype=dtypes)
        except pd.errors.EmptyDataError:
            data = pd.DataFrame(columns=self.columns)
            data = self.cast_dataframe(data)

        if chunksize:
            return return_gen(data)
        else:
            self.__verify_data(data)
            return data


class CsvOutput(object):
    def __init__(
            self, filepath, dtypes, header=True,
            na_rep='NaN', columns=None
    ):
        self.filepath = filepath
        self.header = header
        self.dtypes = dtypes
        self.na_rep = na_rep

        self.columns = columns

        self.__confirm_compression_type_pandas()

        self.sep = ','

    @property
    def yaml_file(self):
        return self.filepath + '.yaml'

    @property
    def header_line(self):
        return self.sep.join(self.columns) + '\n'

    def __confirm_compression_type_pandas(self):
        filepath = self.filepath
        if filepath.endswith('.tmp'):
            filepath = filepath[:-4]

        _, ext = os.path.splitext(filepath)

        if not ext == ".gz":
            raise CsvWriterError("{} is not supported".format(ext))

    def write_yaml(self):
        type_converter = pandas_to_std_types()

        yamldata = {'header': self.header, 'sep': self.sep, 'columns': []}
        for column in self.columns:
            data = {'name': column, 'dtype': type_converter[self.dtypes[column]]}
            yamldata['columns'].append(data)

        with open(self.yaml_file, 'wt') as f:
            yaml.safe_dump(yamldata, f, default_flow_style=False)

    def __get_dtypes_from_df(self, df):
        dtypes = df.dtypes
        dtypes_converter = pandas_to_std_types()
        dtypes = {k: dtypes_converter[str(v)] for k, v in dtypes.items()}
        return dtypes

    def __verify_df(self, df):
        if self.columns:
            assert set(list(df.columns.values)) == set(self.columns)
        else:
            self.columns = df.columns.values

    def __write_df(self, df, header=True, mode='w'):
        self.__verify_df(df)

        df.to_csv(
            self.filepath, sep=self.sep, na_rep=self.na_rep,
            index=False, compression='gzip', mode=mode, header=header
        )

    def __write_df_chunks(self, dfs, header=True):
        for i, df in enumerate(dfs):
            if i == 0 and self.header:
                self.__write_df(df, header=header, mode='w')
            else:
                self.__write_df(df, header=False, mode='a')

    def write_df(self, df, chunks=False):
        if chunks:
            self.__write_df_chunks(df, header=self.header)
        else:
            self.__write_df(df, self.header)

        self.write_yaml()

    def write_header(self, writer):
        header = ','.join(self.columns)
        header = header + '\n'
        writer.write(header)

    def write_data_streams(self, csvfiles):
        assert self.columns
        assert self.dtypes
        with gzip.open(self.filepath, 'wt') as writer:

            if self.header:
                self.write_header(writer)

            for csvfile in csvfiles:
                with gzip.open(csvfile, 'rt') as data_stream:
                    shutil.copyfileobj(
                        data_stream, writer, length=16 * 1024 * 1024
                    )

        self.write_yaml()

    def rewrite_csv(self, csvfile):
        assert self.columns
        assert self.dtypes
        with gzip.open(self.filepath, 'wt') as writer:
            if self.header:
                self.write_header(writer)

            with gzip.open(csvfile, 'rt') as data_stream:
                shutil.copyfileobj(
                    data_stream, writer, length=16 * 1024 * 1024
                )

        self.write_yaml()

    def write_text(self, text):
        assert self.columns
        assert self.dtypes

        with gzip.open(self.filepath, 'wt') as writer:

            if self.header:
                self.write_header(writer)

            for line in text:
                writer.write(line)


def write_metadata(infile):
    csvinput = IrregularCsvInput(infile)

    csvoutput = CsvOutput(
        infile, csvinput.dtypes, header=csvinput.header,
        columns=csvinput.columns
    )
    csvoutput.write_yaml()


def merge_dtypes(dtypes_all):

    if dtypes_all == []:
        raise CsvMergeDtypesEmptyMergeSet("must provide dtypes to merge")

    merged_dtypes = {}

    for dtypes in dtypes_all:
        for k, v in dtypes.items():
            if k in merged_dtypes:
                if merged_dtypes[k] != v:
                    raise DtypesMergeException("dtypes not mergeable")
            else:
                merged_dtypes[k] = v

    return merged_dtypes


def concatenate_csv(inputfiles, output, write_header=True):

    if inputfiles == [] or inputfiles == {}:
        raise CsvConcatException("nothing provided to concat")

    if isinstance(inputfiles, dict):
        inputfiles = inputfiles.values()

    inputs = [CsvInput(infile) for infile in inputfiles]

    dtypes = [csvinput.dtypes for csvinput in inputs]

    #if all types are int
    uq_cols = set(list(itertools.chain(*[d.keys() for d in dtypes])))

    uq_dtypes = set(list(itertools.chain(*[list(d.values()) for d in dtypes])))

    if uq_cols != set(dtypes[0].keys()):
        if uq_dtypes == {"int"}:
            raise CsvConcatNaNIntDtypeException("if dtypes are 'int',"
                                                " all columns must be identical")

    dtypes = merge_dtypes(dtypes)

    headers = [csvinput.header for csvinput in inputs]

    columns = [csvinput.columns for csvinput in inputs]
    columns = list(set([col for column_set in columns for col in column_set]))

    low_memory = True
    if any(headers):
        low_memory = False

    if not all(columns[0] == elem for elem in columns):
        low_memory = False


    if low_memory:
        concatenate_csv_files_quick_lowmem(inputfiles, output, dtypes, columns, write_header=write_header)
    else:
        concatenate_csv_files_pandas(inputfiles, output, dtypes, columns, write_header=write_header)


def concatenate_csv_files_pandas(in_filenames, out_filename, dtypes, columns, write_header=True):

    if isinstance(in_filenames, dict):
        in_filenames = in_filenames.values()

    data = [
        CsvInput(in_filename).read_csv() for in_filename in in_filenames
    ]
    data = pd.concat(data, ignore_index=True)
    csvoutput = CsvOutput(out_filename, dtypes, header=write_header, columns=columns)
    csvoutput.write_df(data)


def concatenate_csv_files_quick_lowmem(inputfiles, output, dtypes, columns, write_header=True):
    csvoutput = CsvOutput(
        output, dtypes, header=write_header, columns=columns
    )
    csvoutput.write_data_streams(inputfiles)


#annotation_dtypes shouldnt be default, if it is None, it breaks
def annotate_csv(infile, annotation_data, outfile, annotation_dtypes, on="cell_id", write_header=True):

    csvinput = CsvInput(infile)
    metrics_df = csvinput.read_csv()

    merge_on = metrics_df[on]

    if not all([on in list(annotation_data.keys()) for on in merge_on]):
        raise CsvAnnotateError("annotation_data must contain all "
                               "elements in infile's annotation col")

    for cell in merge_on:
        col_data = annotation_data[cell]
        for column, value in col_data.items():
            metrics_df.loc[metrics_df[on] == cell, column] = value

    csv_dtypes = csvinput.dtypes

    assert not set(csv_dtypes.keys()).intersection(annotation_dtypes.keys())

    csv_dtypes.update(annotation_dtypes)

    output = CsvOutput(outfile, csv_dtypes, header=write_header)
    output.write_df(metrics_df)


def rewrite_csv_file(filepath, outputfile, write_header=True):
    """
    generate header less csv files
    :param filepath:
    :type filepath:
    :param outputfile:
    :type outputfile:
    """
    csvinput = CsvInput(filepath)

    if csvinput.header:
        if write_header:
            raise Exception('no op')

        df = csvinput.read_csv()

        csvoutput = CsvOutput(
            outputfile, header=write_header, columns=csvinput.columns,
            dtypes=csvinput.dtypes
        )
        csvoutput.write_df(df)

    else:
        if not write_header:
            raise Exception('no op')

        csvoutput = CsvOutput(
            outputfile, header=write_header, columns=csvinput.columns,
            dtypes=csvinput.dtypes
        )
        csvoutput.rewrite_csv(filepath)


def merge_csv(in_filenames, out_filename, how, on, write_header=True):
    if isinstance(in_filenames, dict):
        in_filenames = in_filenames.values()

    data = [CsvInput(infile) for infile in in_filenames]

    dfs = [csvinput.read_csv() for csvinput in data]

    dtypes = [csvinput.dtypes for csvinput in data]

    data = merge_frames(dfs, how, on)

    dtypes = merge_dtypes(dtypes)

    columns = list(data.columns.values)

    csvoutput = CsvOutput(out_filename, dtypes, header=write_header, columns=columns)
    csvoutput.write_df(data)


def _validate_merge_cols(frames, on):
    '''
    make sure frames look good. raise relevant exceptions
    :param frames: list of dfs to merge
    :param on: list of common columns in frames on which to merge
    :return: nothing
    '''
    if on == []:
        raise CsvMergeException("unable to merge if given nothing to merge on")

    #check that columns to be merged have identical values
    standard = frames[0][on]
    for frame in frames:
        if not standard.equals(frame[on]):
            raise CsvMergeColumnMismatchException("columns on which to merge must be identical")

    #check that columns to be merged have same dtypes
    for shared_col in on:
        if len(set([frame[shared_col].dtypes for frame in frames])) != 1:
            raise CsvMergeColumnMismatchException("columns on which to merge must have same dtypes")

    common_cols = set.intersection(*[set(frame.columns) for frame in frames])
    cols_to_check = list(common_cols - set(on))

    for frame1, frame2 in zip(frames[:-1], frames[1:]):
        if not frame1[cols_to_check].equals(frame2[cols_to_check]):
            raise CsvMergeCommonColException("non-merged common cols must be identical")

def merge_frames(frames, how, on):
    """
    annotates input_df using ref_df
    :param frames:
    :type frames:
    :param how:
    :type how:
    :param on:
    :type on:
    :return:
    :rtype:
    """

    if ',' in on:
        on = on.split(',')

    _validate_merge_cols(frames, on)

    if len(frames) == 1:
        return frames[0]

    else:
        left = frames[0]
        right = frames[1]
        cols_to_use = list(right.columns.difference(left.columns))
        cols_to_use += on
        cols_to_use = list(set(cols_to_use))

        merged_frame = pd.merge(
            left, right[cols_to_use], how=how, on=on,
        )
        for i, frame in enumerate(frames[2:]):
            merged_frame = pd.merge(
                merged_frame, frame, how=how, on=on,
            )
        return merged_frame


def write_dataframe_to_csv_and_yaml(df, outfile, dtypes, write_header=True):
    csvoutput = CsvOutput(outfile, dtypes, header=write_header)

    csvoutput.write_df(df)
