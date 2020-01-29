import os
import numpy as np
import pandas as pd
import single_cell.utils.csvutils as csvutils
import random


################################################
###                utilities                ###
################################################
def dfs_match(data, reference, match_levels):
    assert set(data.columns.values) == set(refdata.columns.values)

    for colname in data.columns.values:
        is_exact = match_levels[colname] == "exact":
        
        if is_exact:
            exact_compare_cols(data, refdata, colname)
        else:
            approx_compare_cols(data, reference, colname)

def exact_compare_cols(data, reference, column_name):
    data_index = set(data.index)
    reference_index = set(reference.index)

    assert data_index == reference_index

    index_order = sorted(data_index)
    data = data.reindex(index_order)
    reference = data.reindex(index_order)

    assert data[column_name].equals(reference[column_name])

def approx_compare_cols(data, reference, column_name, eps=0.001):
    data_index = set(data.index)
    reference_index = set(reference.index)

    assert data_index == reference_index

    index_order = sorted(data_index)
    data = data.reindex(index_order)
    reference = data.reindex(index_order)

    diff = data[column_name] - reference[column_name]

    assert np.nanmax(diff.tolist()) < eps

def make_test_df(name, dtypes, tmpdir, length):
    if not name.endswith(".csv"):
        name += ".csv"
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
    df = pd.DataFrame(df_dic, columns = dtypes.keys())
    return filename, csvutils.write_dataframe_to_csv_and_yaml(df, filename, dtypes)

def _rand_float_col(n):
    return [random.uniform(0, 100) for _ in xrange(n)]

def _rand_int_col(n):
    return [random.randint(0, 10) for _ in xrange(n)]

def _rand_bool_col(n):
    return [random.choice([True, False]) for _ in xrange(n)]

def _rand_str_col(n):
    ''.join(random.choices(string.ascii_uppercase + string.digits, k=n))

################################################
###                  tests                   ###
################################################
class TestAnnotateCsv:
    def __init__(self, tmpdir):
        self.tmpdir = tmpdir

    def test_annotate_csv(self):
        '''
        basic sanity check - test annotating normal csv
        '''
        pass 

    def test_annotate_csv_wrong_lengths(self):
        '''
        test annotating csv where annotation_data differs in length from csv
        '''
        pass  

    def test_annotate_csv_dtype_mismatch(self):
        '''
        test annotating csv with inappropriate annotation_dtypes
        '''
        pass     

    def test_annotate_csv_no_write_header(self):
        '''
        test annotating csv without writing header
        '''
        pass           

class TestConcat:
    
    class TestConcatCsv(TestConcat):
        
        def __init__(self, tmpdir):
            self.tmpdir = tmpdir
        
        def test_concat_csv(self):
            '''
            basic sanity check - concat two csvs with same cols
            '''
            pass
        
        def test_concat_csv_no_header(self):
            '''
            concat two csvs without headers
            '''
            pass    
        
        def test_concat_csv_quick_vs_pandas(self):
            '''
            concat two csvs using both methods
            '''
            pass
        
    class TestConcatCsvFilesPandas(TestConcat):
        def __init__(self, tmpdir):
            self.tmpdir = tmpdir

        def test_concat_csv(self):
            '''
            basic sanity check - concat 1 csvs with same cols
            '''
            pass   

        def test_concat_csv_nothing_to_concat(self):
            '''
            provide nothing to concat
            '''
            pass    
        
        def test_concat_csv_different_cols(self):
            '''
            concat two dataframes with different columns
            '''
            pass    
        
        def test_concat_csv_one_file_to_concat(self):
            '''
            provide just 1 file to concat
            '''
            pass   
        
        def test_concat_csv_different_dtypes(self):
            '''
            concat two dataframes same colnames with different dtypes
            '''
            pass    

        def test_concat_csv_with_nans(self):
            '''
            concat two csvs with NaNs
            '''
            pass    
    
    class TestConcatCsvFilesQuickLowmem(TestConcat):
        def __init__(self, tmpdir):
            self.tmpdir = tmpdir
        
        def test_quick_concat(self):
            '''
            sanity check - two easy csvs
            '''
            pass    

        def test_quick_concat_scale(self):
            '''
            concat two large dataframes
            '''
            pass

        def test_quick_concat_csv_nothing_to_concat(self):
            '''
            provide nothing to concat
            '''
            pass    

        def test_quick_concat_csv_different_cols(self):
            '''
            concat two dataframes with different columns
            '''
            pass    

        def test_quick_concat_csv_one_file_to_concat(self):
            '''
            provide just 1 file to concat
            '''
            pass   

        def test_quick_concat_csv_different_dtypes(self):
            '''
            concat two dataframes same colnames with different dtypes
            '''
            pass    

        def test_quick_concat_csv_with_nans(self):
            '''
            concat two csvs with NaNs
            '''
            pass    

class TestWrite:

    class TestWriteMetadata(TestWrite):
        
        def __init__(self, tmpdir):
            self.tmpdir = tmpdir
        
        def test_write_metadata(self):
            '''
            basic sanity check - test writing metadata
            '''
            pass    

        def test_write_metadata_emty_input(self):
            '''
            test writing metadata with empty csv
            '''
            pass    

        def test_write_metadata_no_header(self):
            '''
            test writing metadata with csv without header
            '''
            pass    

        def test_write_metadata_with_NaNs(self):
            '''
            basic sanity check - test writing metadata
            '''
            pass          
        
    class TestRewriteCsvFile(TestWrite):
        
        def __init__(self, tmpdir):
            self.tmpdir = tmpdir

        def test_rewrite_csv_file(self):
            '''
            basic sanity check - test rewriting  csv file
            '''
            pass    

        def test_rewrite_csv_file_no_write_header(self):
            '''
            test rewriting  csv file without writing header
            '''
            pass    

        def test_rewrite_csv_file_no_header_write_anyway(self):
            '''
            test rewriting  csv file without header but 
            force function to write out anyway
            '''
            pass      
    
        def test_rewrite_csv_file_empty_csv(self):
            '''
            test rewriting empty csv file
            '''
            pass    

    class WriteDataframeToCsvAndYaml(TestWrite):
        def __init__(self, tmpdir):
            self.tmpdir = tmpdir
        
        def test_write_to_csv_yaml(self):
            '''
            basic sanity check - write normal df
            '''
            pass   

        def test_write_to_csv_yaml_no_header(self):
            '''
            write single df without header
            '''
            pass   

        def test_write_to_csv_yaml_empty(self):
            '''
            write empty df
            '''
            pass   

class TestMerge:
    
    class TestMergeFrames(TestMerge):
        def __init__(self, tmpdir):
            self.tmpdir = tmpdir

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
        
        def test_merge_multiple_cols(self, n_cols):
            """
            test merging of 2 dfs on multiple columns with right merge
            """
            pass
        
        def test_merge_frames_multiple_frames(self, n_frames):
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
    
    class TestMergeDtypes(TestMerge):
        
        def __init__(self, tmpdir):
            self.tmpdir = tmpdir
            
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
    
    class TestMergeCsv(TestMerge):

        def __init__(self, tmpdir):
            self.tmpdir = tmpdir

        def test_merge_csv(self):
            dtypes1 = {v: str(int) for v in 'ABCD'}
            name1, df1 = make_test_df('first.csv', dtypes1, self.tmpdir, 100):

            dtypes2 = {v: str(int) for v in 'AEFGH'}
            name1, df1 = make_test_df('second.csv', dtypes2, self.tmpdir, 100):

            merged = os.path.join(self.tmpdir, 'merged.csv')

            df2['A'] = df1['A']

            csvutils.write_dataframe_to_csv_and_yaml(df1, csv_1, dtypes)
            csvutils.write_dataframe_to_csv_and_yaml(df2, csv_2, dtypes)

            csvutils.merge_csv([csv_1, csv_2], merged, 'outer', ['A'])

        def test_merge_csv_two_cols(self):
            """
            :param tmpdir:
            :type tmpdir:
            :return:
            :rtype:
            """
            pass

        def test_merge_csv_gzipped(self):
            pass

        def test_merge_csv_input_format(self):
            pass

        def test_merge_csv_no_header(self):
            pass


# if colvals:
#             for col, vals in colvals.items():
#                 assert vals.dtype == dtypes[col]
#                 df_dict[col] = vals
#         random_cols = list(set(dtypes.keys()) - set(colvals.keys()))
#         if not any random_cols:
#             return