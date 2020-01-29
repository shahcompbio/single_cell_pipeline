import os

import numpy as np
import pandas as pd
import single_cell.utils.csvutils as csvutils



class MergeCsvTests:

    def test_merge_csv(tmpdir):
        csv_1 = os.path.join(tmpdir, 'first.csv')
        csv_2 = os.path.join(tmpdir, 'second.csv')
        merged = os.path.join(tmpdir, 'merged.csv')

        df1 = pd.DataFrame(np.random.randint(0, 100, size=(100, 4)), columns=list('ABCD'))
        df2 = pd.DataFrame(np.random.randint(0, 100, size=(100, 4)), columns=list('EFGH'))

        df2['A'] = df1['A']

        dtypes = {v: str(int) for v in 'ABCD'}
        csvutils.write_dataframe_to_csv_and_yaml(df1, csv_1, dtypes)

        dtypes = {v: str(int) for v in 'AEFGH'}
        csvutils.write_dataframe_to_csv_and_yaml(df2, csv_2, dtypes)

        csvutils.merge_csv([csv_1, csv_2], merged, 'outer', ['A'])


    def test_merge_csv_two_cols(tmpdir):
        """
        :param tmpdir:
        :type tmpdir:
        :return:
        :rtype:
        """
        pass
