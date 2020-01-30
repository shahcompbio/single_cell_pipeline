import os
import numpy as np
import pandas as pd
import single_cell.utils.csvutils as csvutils
import random
import string
import pytest

###############################################
#                fixtures                     #
################################################


@pytest.fixture
def n_frames():
    return random.randint(2, 10)

###############################################
#                utilities                    #
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
            df_dict[col] = "".join(_enumerated_list(length))
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
    return range(n)
    #return [random.randint(0, 10) for _ in range(n)]


def _rand_bool_col(n):
    return [random.choice([True, False]) for _ in range(n)]


def _enumerated_list(n, must_have="", count=0):
    s = random.choices(string.ascii_uppercase, k=n)
    if count:
        s = [l+str(count) for l in s]
    s.append(must_have)
    return s



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


class TestMergeHelpers:

    def make_mergeable_test_dfs(self, tmpdir, n_dfs, dtypes,
                                shared, length, write=False,
                                headers=True, get_expected=False,
                                merge_args=None):

        dfs = []
        names = []

        for dtype_set in dtypes:
            assert set(shared).issubset(set(dtype_set.keys()))

        for i in range(n_dfs):
            base_name = str(i) + ".csv.gz"
            name, df = make_test_df(base_name, dtypes[i], tmpdir, length)
            names.append(name)
            dfs.append(df)

        shared_values = dfs[0][shared] #arbitary choose first df

        for df in dfs:
            df.update(shared_values)

        expected_output = None
        if write:
            for i in range(n_dfs):
                csvutils.write_dataframe_to_csv_and_yaml(dfs[i],
                                                         names[i],
                                                         dtypes[i],
                                                         headers)
        print ("i")
        if get_expected and merge_args:
            print (dfs[0], dfs[1], merge_args["suffixes"])
            expected_output = dfs[0].merge(dfs[1], how=merge_args["how"],
                                           on=merge_args["on"],
                                           suffixes=merge_args["suffixes"])
        print( "a")
        return dfs, names, expected_output


class TestMergeFrames(TestMergeHelpers):

    def test_merge_frames(self):
        """
        basic sanity check - test merging of two frames on 1 col
        """

        dtypes1 = {v: "int" for v in 'ABCD'}
        dtypes2 = {v: "int" for v in 'AEFGH'}

        merge_args = {"how": 'outer', "on": ['A'], "suffixes": ["", ""]}

        dfs, names, ref = self.make_mergeable_test_dfs("", 2,
                                                       [dtypes1, dtypes2], ["A"],
                                                       100, write=False,
                                                       get_expected=True,
                                                       merge_args=merge_args)

        merged = csvutils.merge_frames(dfs, how=merge_args["how"],
                                       on=merge_args["on"],
                                       suffixes=merge_args["suffixes"])

        assert dfs_exact_match(ref, merged)

    def test_merge_frames_inner(self):
        """
        test merging of 2 dfs on 1 col with inner merge
        """

        dtypes1 = {v: "int" for v in 'ABCD'}
        dtypes2 = {v: "int" for v in 'AEFGH'}

        merge_args = {"how": 'inner', "on": ['A'], "suffixes": ["", ""]}

        dfs, names, ref = self.make_mergeable_test_dfs("", 2,
                                                       [dtypes1, dtypes2], ["A"],
                                                       100, write=False,
                                                       get_expected=True,
                                                       merge_args=merge_args)

        merged = csvutils.merge_frames(dfs, how=merge_args["how"],
                                       on=merge_args["on"],
                                       suffixes=merge_args["suffixes"])

        assert dfs_exact_match(ref, merged)

    def test_merge_frames_left(self):
        """
        test merging of 2 dfs on 1 col with left merge
        """

        dtypes1 = {v: "int" for v in 'ABCD'}
        dtypes2 = {v: "int" for v in 'AEFGH'}

        merge_args = {"how": 'left', "on": ['A'], "suffixes": ["", ""]}

        dfs, names, ref = self.make_mergeable_test_dfs("", 2,
                                                       [dtypes1, dtypes2], ["A"],
                                                       100, write=False,
                                                       get_expected=True,
                                                       merge_args=merge_args)

        merged = csvutils.merge_frames(dfs, how=merge_args["how"],
                                       on=merge_args["on"],
                                       suffixes=merge_args["suffixes"])

        assert dfs_exact_match(ref, merged)

    def test_merge_frames_right(self):
        """
        test merging of 2 dfs on 1 col with right merge
        """

        dtypes1 = {v: "int" for v in 'ABCD'}
        dtypes2 = {v: "int" for v in 'AEFGH'}

        merge_args = {"how": 'right', "on": ['A'], "suffixes": ["", ""]}

        dfs, names, ref = self.make_mergeable_test_dfs("", 2,
                                                       [dtypes1, dtypes2], ["A"],
                                                       100, write=False,
                                                       get_expected=True,
                                                       merge_args=merge_args)

        merged = csvutils.merge_frames(dfs, how=merge_args["how"],
                                       on=merge_args["on"],
                                       suffixes=merge_args["suffixes"])

        assert dfs_exact_match(ref, merged)

    def test_merge_multiple_cols(self):
        """
        test merging of 2 dfs on multiple columns with right merge
        """
        pass

    def test_merge_frames_multiple_frames(self, n_frames):
        """
        test merging of n_frames on 1 col with right merge
        :param n_frames: number of dataframes to merge
        """

        dtypes = []
        for i in range(n_frames):
            dtypes.append({v: "int" for v in _enumerated_list(4,
                                                              must_have="A",
                                                              count=i)})

        print(dtypes)
        merge_args = {"how": 'right', "on": ['A'], "suffixes": [""]*n_frames}

        print("d")
        dfs, names, ref = self.make_mergeable_test_dfs("", n_frames, dtypes,
                                                       ["A"], 100, write=False,
                                                       get_expected=True,
                                                       merge_args=merge_args)
        print ("here")
        print (dfs)
        merged = csvutils.merge_frames(dfs, how=merge_args["how"],
                                       on=merge_args["on"],
                                       suffixes=merge_args["suffixes"])

        print("f")
        assert dfs_exact_match(ref, merged)

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


class TestMergeDtypes(TestMergeHelpers):

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


class TestMergeCsv(TestMergeHelpers):

    def test_merge_csv(self, tmpdir):

        dtypes1 = {v: "int" for v in 'ABCD'}
        dtypes2 = {v: "int" for v in 'AEFGH'}

        merged = os.path.join(tmpdir, 'merged.csv.gz')
        merge_args = {"how":'outer', "on":['A'], "suffixes":["",""]}

        dfs, names, ref = self.make_mergeable_test_dfs(tmpdir, 2,
                                                  [dtypes1, dtypes2], ["A"],
                                                  100, write=True,
                                                  get_expected=True,
                                                  merge_args=merge_args)

        csvutils.merge_csv(names, merged, how=merge_args["how"],
                           on=merge_args["on"],
                           suffixes=merge_args["suffixes"])

        assert os.path.exists(merged)

        assert dfs_exact_match(ref, merged)


    def test_merge_csv_no_header(self, tmpdir):
        """
        merge csvs with no headers
        """

        dtypes1 = {v: "int" for v in 'ABCD'}
        dtypes2 = {v: "int" for v in 'AEFGH'}

        merged = os.path.join(tmpdir, 'merged.csv.gz')
        merge_args = {"how": 'outer', "on": ['A'], "suffixes": ["", ""]}

        dfs, names, ref = self.make_mergeable_test_dfs(tmpdir, 2,
                                                       [dtypes1, dtypes2], ["A"],
                                                       100, headers=False,
                                                       write=True,
                                                       get_expected=True,
                                                       merge_args=merge_args)

        csvutils.merge_csv(names, merged, how=merge_args["how"],
                           on=merge_args["on"],
                           suffixes=merge_args["suffixes"])

        assert os.path.exists(merged)

        assert dfs_exact_match(ref, merged)

