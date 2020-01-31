import os
import numpy as np
import pandas as pd
import single_cell.utils.csvutils as csvutils
import random
import string
import pytest
import itertools

###############################################
#                fixtures                     #
################################################

@pytest.fixture
def n_dtypes():
    return random.randint(2, 10)

@pytest.fixture
def n_frames():
    return random.randint(2, 10)

@pytest.fixture
def n_rows():
    return 10
###############################################
#                utilities                    #
################################################


def _raises_correct_error(function, *args,
                          expected_error=csvutils.CsvParseError,
                          **kwargs):

    raised = False
    try:
        function(*args, **kwargs)
    except Exception as e:
        if type(e) == expected_error:
            raised = True
        else:
            print("raised wrong error: raised: {}, expected: {}"
                   .format(type(e), expected_error))
    finally:
        return raised

def dfs_exact_match(data, reference):

    if isinstance(data, str):
        data = csvutils.CsvInput(data).read_csv()
    if isinstance(reference, str):
        reference = csvutils.CsvInput(reference).read_csv()
    if set(data.columns) != set(reference.columns):
        return False
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
            df_dict[col] = _str_list(length)
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


def _rand_int_col(n, scale=1):
    return list(range(n)) * scale
    #return [random.randint(0, 10) for _ in range(n)]


def _rand_bool_col(n):
    return [random.choice([True, False]) for _ in range(n)]


def _str_list(n, must_have="", count=0):
    s = random.sample(string.ascii_uppercase, k=n)
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

    #these take length so as to not tap into the fixture above
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

        if get_expected and merge_args:
            expected_output = dfs[0].merge(dfs[1], how=merge_args["how"],
                                           on=merge_args["on"],
                                           suffixes=merge_args["suffixes"])

        return dfs, names, expected_output

    def merge_directional_test(self, length, direction):

        dtypes1 = {v: "int" for v in 'ABCD'}
        dtypes2 = {v: "int" for v in 'AEFGH'}

        merge_args = {"how": direction, "on": ['A'], "suffixes": ["", ""]}

        dfs, names, ref = self.make_mergeable_test_dfs("", 2,
                                                       [dtypes1, dtypes2], ["A"],
                                                       length, write=False,
                                                       get_expected=True,
                                                       merge_args=merge_args)

        merged = csvutils.merge_frames(dfs, how=merge_args["how"],
                                       on=merge_args["on"])

        assert dfs_exact_match(ref, merged)

class TestMergeFrames(TestMergeHelpers):

    def test_merge_frames(self, n_rows):
        """
        basic sanity check - test merging of two frames on 1 col
        """

        self.merge_directional_test(n_rows, "outer")

    def test_merge_frames_inner(self, n_rows):
        """
        test merging of 2 dfs on 1 col with inner merge
        """
        self.merge_directional_test(n_rows, "inner")

    def test_merge_frames_left(self, n_rows):
        """
        test merging of 2 dfs on 1 col with left merge
        """

        self.merge_directional_test(n_rows, "left")

    def test_merge_frames_right(self, n_rows):
        """
        test merging of 2 dfs on 1 col with right merge
        """

        self.merge_directional_test(n_rows, "right")

    def test_merge_frames_multiple_cols(self, n_rows):
        """
        test merging of 2 dfs on multiple columns with right merge
        """

        dtypes1 = {v: "int" for v in 'ABCD'}
        dtypes2 = {v: "int" for v in 'ABEFGH'}

        merge_args = {"how": "outer", "on": ['A', 'B'], "suffixes": ["", ""]}

        dfs, names, ref = self.make_mergeable_test_dfs("", 2,
                                                       [dtypes1, dtypes2], merge_args["on"],
                                                       n_rows, write=False,
                                                       get_expected=True,
                                                       merge_args=merge_args)

        merged = csvutils.merge_frames(dfs, how=merge_args["how"],
                                       on=merge_args["on"])
        assert dfs_exact_match(ref, merged)

    def test_merge_frames_one_frame(self, n_rows):
        '''
        provide just one df
        :param n_rows: number of rows in simulated df
        :return: assertion
        '''
        dtypes1 = {v: "int" for v in 'ACD'}

        merge_args = {"how": "outer", "on": ["A"], "suffixes": ["", ""]}

        dfs, name, _ = self.make_mergeable_test_dfs("", 1, [dtypes1],
                                                   merge_args["on"],
                                                   n_rows, write=False,
                                                   merge_args=merge_args) #returns list of 1 df

        merged = csvutils.merge_frames(dfs, how=merge_args["how"],
                                       on=merge_args["on"])

        assert dfs_exact_match(dfs[0], merged)

    def test_merge_frames_multiple_frames(self, n_frames, n_rows):
        """
        test merging of n_frames on 1 col with right merge
        :param n_frames: number of dataframes to merge
        """

        dtypes = []
        n_cols_per_frame = 2
        for i in range(n_frames):
            dtypes.append({v: "int" for v in _str_list(n_cols_per_frame,
                                                              must_have="merge_on",
                                                              count=i)})

        merge_args = {"how": 'right', "on": ['merge_on'], "suffixes": [""]*n_frames}

        dfs, names, _ = self.make_mergeable_test_dfs("", n_frames, dtypes,
                                                       ["merge_on"], n_rows, write=False,
                                                       get_expected=False,
                                                       merge_args=merge_args)

        merged = csvutils.merge_frames(dfs, how=merge_args["how"],
                                       on=merge_args["on"])

        #just naively make sure it has right # of cols
        assert merged.shape[0] == n_rows \
               and merged.shape[1] == (n_cols_per_frame * n_frames) + 1

    def test_merge_frames_nothing_to_merge(self, n_rows):
        """
        test merging of 2 dfs on 0 col
        """

        dtypes1 = {v: "int" for v in 'CD'}
        dtypes2 = {v: "int" for v in 'EFGH'}

        merge_args = {"how": "outer", "on": [], "suffixes": ["", ""]}

        dfs, names, ref = self.make_mergeable_test_dfs("", 2, [dtypes1, dtypes2],
                                                       merge_args["on"],
                                                       n_rows,
                                                       write=False, get_expected=False,
                                                       merge_args=merge_args)

        assert _raises_correct_error(csvutils.merge_frames, dfs,
                                     how=merge_args["how"],
                                     on=merge_args["on"],
                                     expected_error=csvutils.CsvMergeException)

    def test_merge_frames_cols_to_merge_have_different_dtypes(self, n_rows):
        """
        test merging of 2 dfs on 1 col. Each dataframe has different dtype for col
        """

        dtypes1 = {v: "int" for v in 'ABCD'}
        dtypes2 = {v: "str" for v in 'AEFGH'}

        merge_args = {"how": "outer", "on": ['A'], "suffixes": ["", ""]}

        dfs, names, _ = self.make_mergeable_test_dfs("", 2,
                                                     [dtypes1, dtypes2], ["A"],
                                                     n_rows, write=False,
                                                     get_expected=False,
                                                     merge_args=merge_args)

        assert _raises_correct_error(csvutils.merge_frames, dfs,
                                     how=merge_args["how"],
                                     on=merge_args["on"],
                                     expected_error=csvutils.CsvMergeColumnMismatchException)

    def test_merge_frames_create_nans(self, n_rows):
        """
        test merging of 2 dfs on 1 col to create NaNs
        """

        dtypes1 = {v: "int" for v in 'ACD'}
        dtypes2 = {v: "int" for v in 'ACD'}

        merge_args = {"how": "outer", "on": ["A"], "suffixes": ["", ""]}

        dfs, name, _ = self.make_mergeable_test_dfs("", 2, [dtypes1, dtypes2],
                                                    [], n_rows, write=False,
                                                    merge_args=merge_args)
        on = merge_args["on"]
        dfs[0].iloc[2, dfs[0].columns.get_loc(on[0])] = 10

        assert _raises_correct_error(csvutils.merge_frames, dfs,
                                     how=merge_args["how"],
                                     on=merge_args["on"],
                                     expected_error=csvutils.CsvMergeColumnMismatchException)

    def test_merge_frames_with_nans(self, n_rows):
        """
        test merging of 2 dfs on 1 col which contains NaNs in each
        """

        dtypes1 = {v: "float" for v in 'ACD'}
        dtypes2 = {v: "float" for v in 'ACD'}

        merge_args = {"how": "outer", "on": ["A"], "suffixes": ["", ""]}

        dfs, name, _ = self.make_mergeable_test_dfs("", 2, [dtypes1, dtypes2],
                                                    [], n_rows, write=False,
                                                    merge_args=merge_args)
        on = merge_args["on"]
        dfs[0].iloc[2, dfs[0].columns.get_loc(on[0])] = np.NaN
        dfs[1].iloc[2, dfs[1].columns.get_loc(on[0])] = np.NaN

        assert _raises_correct_error(csvutils.merge_frames, dfs,
                                     how=merge_args["how"],
                                     on=merge_args["on"],
                                     expected_error=csvutils.CsvMergeColumnMismatchException)


class TestMergeDtypes(TestMergeHelpers):

    def test_merge_dtypes(self):
        """
        basic sanity check - test merging of two dtype dicts
        """
        dtypes1 = {v: "float" for v in 'ACD'}
        dtypes2 = {v: "float" for v in 'ACDEF'}
        ref = {v: "float" for v in
               set(dtypes1.keys()).union(set(dtypes2.keys()))}

        merged_dtypes = csvutils.merge_dtypes([dtypes1, dtypes2])

        assert ref == merged_dtypes

    def test_merge_dtypes_none_given(self):
        """
        test merging of empty list of dtypes
        """

        assert _raises_correct_error(csvutils.merge_dtypes, [],
                                     expected_error=csvutils.CsvMergeDtypesEmptyMergeSet)

    def test_merge_dtypes_one_given(self):
        """
        test merging of list of 1 dtype dict
        """
        dtypes1 = {v: "float" for v in 'ACD'}
        ref = dtypes1

        merged_dtypes = csvutils.merge_dtypes([dtypes1])

        assert ref == merged_dtypes

    def test_merge_dtypes_multiple_given(self, n_dtypes):
        """
        test merging of n_dtypes dtype dicts
        :param n_dtypes: number of dtypes to merge
        """
        dtypes = [{v: "int" for v in "".join(_str_list(3, "A"))}
                  for _ in range(n_dtypes)]

        #taken from https://stackoverflow.com/questions/9819602/union-of-dict-objects-in-python
        ref = dict(itertools.chain.from_iterable(dct.items()
                                                 for dct in dtypes))

        merged_dtypes = csvutils.merge_dtypes(dtypes)

        assert ref == merged_dtypes

class TestMergeCsv(TestMergeHelpers):

    def test_merge_csv(self, tmpdir, n_rows):

        dtypes1 = {v: "int" for v in 'ABCD'}
        dtypes2 = {v: "int" for v in 'AEFGH'}

        merged = os.path.join(tmpdir, 'merged.csv.gz')
        merge_args = {"how":'outer', "on":['A'], "suffixes":["",""]}

        dfs, names, ref = self.make_mergeable_test_dfs(tmpdir, 2,
                                                  [dtypes1, dtypes2], ["A"],
                                                  n_rows, write=True,
                                                  get_expected=True,
                                                  merge_args=merge_args)

        csvutils.merge_csv(names, merged, how=merge_args["how"],
                           on=merge_args["on"])

        assert os.path.exists(merged)

        assert dfs_exact_match(ref, merged)


    def test_merge_csv_no_header(self, tmpdir, n_rows):
        """
        merge csvs with no headers
        """

        dtypes1 = {v: "int" for v in 'ABCD'}
        dtypes2 = {v: "int" for v in 'AEFGH'}

        merged = os.path.join(tmpdir, 'merged.csv.gz')
        merge_args = {"how": 'outer', "on": ['A'], "suffixes": ["", ""]}

        dfs, names, ref = self.make_mergeable_test_dfs(tmpdir, 2,
                                                       [dtypes1, dtypes2], ["A"],
                                                       n_rows, headers=False,
                                                       write=True,
                                                       get_expected=True,
                                                       merge_args=merge_args)

        csvutils.merge_csv(names, merged, how=merge_args["how"],
                           on=merge_args["on"])

        assert os.path.exists(merged)

        assert dfs_exact_match(ref, merged)

