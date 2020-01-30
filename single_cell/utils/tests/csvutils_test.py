import os
import numpy as np
import pandas as pd
import single_cell.utils.csvutils as csvutils
import random
import string
import pytest

################################################
#                utilities                     #
################################################


def dfs_exact_match(data, reference):
    if isinstance(data, str):
        data = csvutils.CsvInput(data).read_csv()
    if isinstance(reference, str):
        reference = csvutils.CsvInput(reference).read_csv()
    print ("using pandas", data,"using csvutils", reference)
    return all([data[col].equals(reference[col]) for col in reference.columns])

def dfs_match(data, reference, match_levels):
    assert set(data.columns.values) == set(reference.columns.values)


    for col_name in data.columns.values:
        is_exact = match_levels[col_name] == "exact"

        if is_exact:
            exact_compare_cols(data, reference, col_name)
        else:
            approx_compare_cols(data, reference, col_name)


def reset_indexes(data, reference):
    data_index = set(data.index)

    reference_index = set(reference.index)

    assert data_index == reference_index

    index_order = sorted(data_index)
    data = data.reindex(index_order)
    reference = data.reindex(index_order)

    return data, reference


def exact_compare_cols(data, reference, column_name):
    data, reference = reset_indexes(data, reference)

    assert data[column_name].equals(reference[column_name])


def approx_compare_cols(data, reference, column_name, eps=0.001):
    data, reference = reset_indexes(data, reference)

    diff = data[column_name] - reference[column_name]

    assert np.nanmax(diff.tolist()) < eps


def make_test_df(name, dtypes, tmpdir, length, write = False):
    filename = os.path.join(tmpdir, name)
    df_dict = {}
    for col, dtype in dtypes.items():
        if dtype == "str":
            df_dict[col] = _rand_str_col(length)
        if dtype == "int":
            df_dict[col] = _rand_int_col(length)
        if dtype == "float":
            df_dict[col] = _rand_float_col(length)
        if dtype == "bool":
            df_dict[col] = _rand_bool_col(length)
    df = pd.DataFrame(df_dict, columns=dtypes.keys())
    if write:
        csvutils.write_dataframe_to_csv_and_yaml(df, filename, dtypes)

    return filename, df


def _rand_float_col(n):
    return [random.uniform(0, 100) for _ in range(n)]


def _rand_int_col(n):
    return [random.randint(0, 10) for _ in range(n)]


def _rand_bool_col(n):
    return [random.choice([True, False]) for _ in range(n)]


def _rand_str_col(n):
    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=n))


################################################
#                  tests                       #
################################################

class TestAnnotateCsv:

    def test_annotate_csv(self):
        """
        basic sanity check - test annotating normal csv
        """
        pass

    def test_annotate_csv_wrong_lengths(self):
        """
        test annotating csv where annotation_data differs in length from csv
        """
        pass

    def test_annotate_csv_dtype_mismatch(self):
        """
        test annotating csv with inappropriate annotation_dtypes
        """
        pass

    def test_annotate_csv_no_write_header(self):
        """
        test annotating csv without writing header
        """
        pass


class TestConcatCsv:

    def test_concat_csv(self):
        """
        basic sanity check - concat two csvs with same cols
        """
        pass

    def test_concat_csv_no_header(self):
        """
        concat two csvs without headers
        """
        pass

    def test_concat_csv_quick_vs_pandas(self):
        """
        concat two csvs using both methods
        """
        pass


class TestConcatCsvFilesPandas:

    def test_concat_csv(self):
        """
        basic sanity check - concat 1 csvs with same cols
        """
        pass

    def test_concat_csv_nothing_to_concat(self):
        """
        provide nothing to concat
        """
        pass

    def test_concat_csv_different_cols(self):
        """
        concat two dataframes with different columns
        """
        pass

    def test_concat_csv_one_file_to_concat(self):
        """
        provide just 1 file to concat
        """
        pass

    def test_concat_csv_different_dtypes(self):
        """
        concat two dataframes same colnames with different dtypes
        """
        pass

    def test_concat_csv_with_nans(self):
        """
        concat two csvs with NaNs
        """
        pass


class TestConcatCsvFilesQuickLowMem:

    def test_quick_concat(self):
        """
        sanity check - two easy csvs
        """
        pass

    def test_quick_concat_scale(self):
        """
        concat two large dataframes
        """
        pass

    def test_quick_concat_csv_nothing_to_concat(self):
        """
        provide nothing to concat
        """
        pass

    def test_quick_concat_csv_different_cols(self):
        """
        concat two dataframes with different columns
        """
        pass

    def test_quick_concat_csv_one_file_to_concat(self):
        """
        provide just 1 file to concat
        """
        pass

    def test_quick_concat_csv_different_dtypes(self):
        """
        concat two dataframes same colnames with different dtypes
        """
        pass

    def test_quick_concat_csv_with_nans(self):
        """
        concat two csvs with NaNs
        """
        pass


class TestWriteMetadata:

    def test_write_metadata(self):
        """
        basic sanity check - test writing metadata
        """
        pass

    def test_write_metadata_emty_input(self):
        """
        test writing metadata with empty csv
        """
        pass

    def test_write_metadata_no_header(self):
        """
        test writing metadata with csv without header
        """
        pass

    def test_write_metadata_with_NaNs(self):
        """
        basic sanity check - test writing metadata
        """
        pass


class TestRewriteCsvFile:

    def test_rewrite_csv_file(self):
        """
        basic sanity check - test rewriting  csv file
        """
        pass

    def test_rewrite_csv_file_no_write_header(self):
        """
        test rewriting  csv file without writing header
        """
        pass

    def test_rewrite_csv_file_no_header_write_anyway(self):
        """
        test rewriting  csv file without header but
        force function to write out anyway
        """
        pass

    def test_rewrite_csv_file_empty_csv(self):
        """
        test rewriting empty csv file
        """
        pass


class WriteDataFrameToCsvAndYaml:

    def test_write_to_csv_yaml(self):
        """
        basic sanity check - write normal df
        """
        pass

    def test_write_to_csv_yaml_no_header(self):
        """
        write single df without header
        """
        pass

    def test_write_to_csv_yaml_empty(self):
        """
        write empty df
        """
        pass


class TestMergeFrames:

    def test_merge_frames(self):
        """
        basic sanity check - test merging of two frames on 1 col
        """
        pass

    def test_merge_frames_inner(self):
        """
        test merging of 2 dfs on 1 col with inner merge
        """
        pass

    def test_merge_frames_outer(self):
        """
        test merging of 2 dfs on 1 col with outer merge
        """
        pass

    def test_merge_frames_left(self):
        """
        test merging of 2 dfs on 1 col with left merge
        """
        pass

    def test_merge_frames_right(self):
        """
        test merging of 2 dfs on 1 col with right merge
        """
        pass

    def test_merge_frames_no_frames(self):
        """
        test merging of 0 dfs
        """
        pass

    def test_merge_frames_one_frame(self):
        """
        test merging of 1 df on 1 col with right merge
        """
        pass

    def test_merge_multiple_cols(self):
        """
        test merging of 2 dfs on multiple columns with right merge
        """
        pass

    def test_merge_frames_multiple_frames(self):
        """
        test merging of n_frames on 1 col with right merge
        :param n_frames: number of dataframes to merge
        """
        pass

    def test_merge_frames_nothing_to_merge(self):
        """
        test merging of 2 dfs on 0 col
        """
        pass

    def test_merge_frames_cols_to_merge_have_different_dtypes(self):
        """
        test merging of 2 dfs on 1 col. Each dataframe has different dtype for col
        """
        pass

    def test_merge_frames_create_nans(self):
        """
        test merging of 2 dfs on 1 col to create NaNs
        """
        pass

    def test_merge_frames_with_nans(self):
        """
        test merging of 2 dfs on 1 col which contains NaNs in each
        """
        pass


class TestMergeDtypes:

    def test_merge_dtypes(self):
        """
        basic sanity check - test merging of two dtype dicts
        """
        pass

    def test_merge_dtypes_none_given(self):
        """
        test merging of empty list of dtypes
        """
        pass

    def test_merge_dtypes_one_give(self):
        """
        test merging of list of 1 dtype dict
        """
        pass

    def test_merge_dtypes_multiple_given(n_dtypes):
        """
        test merging of n_dtypes dtype dicts
        :param n_dtypes: number of dtypes to merge
        """
        pass


class TestMergeCsv:

    def test_merge_csv(self, tmpdir):

        dtypes1 = {v: "int" for v in 'ABCD'}
        name1, df1 = make_test_df('first.csv.gz', dtypes1, tmpdir, 100)

        dtypes2 = {v: "int" for v in 'AEFGH'}
        name2, df2 = make_test_df('second.csv.gz', dtypes2, tmpdir, 100)

        merged = os.path.join(tmpdir, 'merged.csv.gz')

        df2['A'] = df1['A']
        assert df2['A'].equals(df1['A'])
        print ("df1", df1, "df2", df2)

        csvutils.write_dataframe_to_csv_and_yaml(df1, name1, dtypes1)
        csvutils.write_dataframe_to_csv_and_yaml(df2, name2, dtypes2)

        csvutils.merge_csv([name1, name2], merged, 'outer', ['A'], suffixes=["",""])

        assert os.path.exists(merged)

        expected_output = df1.merge(df2, how='outer', on=['A'], suffixes=["",""])

        assert dfs_exact_match(expected_output, merged)

        #errror


    def test_merge_csv_multiple_cols(self):
        """

        """
        pass

    def test_merge_csv_input_format(self):
        """
        """
        pass

    def test_merge_csv_no_header(self):
        """
        """
        pass

