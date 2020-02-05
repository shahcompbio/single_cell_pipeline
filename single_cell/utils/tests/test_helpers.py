import string
import single_cell.utils.csvutils as csvutils
import pandas as pd
import random
import os

class TestInputs:
    #def make_concatable_test_dfs(self, tmpdir, dtypes, concat_args, length):

    def write_dfs(self, tmpdir, dfs, dtypes, write_heads=True):
        n_dfs = len(dfs)
        names = [os.path.join(tmpdir, str(i) + ".csv.gz") for i in range(n_dfs)]

        assert len({n_dfs, len(dtypes)}) == 1

        for i in range(n_dfs):
            csvutils.write_dataframe_to_csv_and_yaml(dfs[i], names[i],
                                                     dtypes[i], write_heads)
        return names

    def make_test_df(self, dtypes, length):

        df_dict = {}
        for col, dtype in dtypes.items():
            df_dict[col] = self.simulate_col(dtype, length)

        df = pd.DataFrame(df_dict, columns=dtypes.keys())
        df = df.astype(dtypes)

        return df

    def make_test_dfs(self, dtypes, length):

        n_dfs = len(dtypes)
        dfs = [self.make_test_df(dtypes[i], length) for i in range(n_dfs)]

        return dfs

    def simulate_col(self, dtype, length):
        if dtype == "str":
            return self._str_list(length)
        if dtype == "int":
            return self._rand_int_col(length)
        if dtype == "float":
            return self._rand_float_col(length)
        if dtype == "bool":
            return self._rand_bool_col(length)

    def _rand_float_col(self, n):
        return [random.uniform(0, 100) for _ in range(n)]

    def _rand_int_col(self, n, scale=1):
        return list(range(n)) * scale
        #return [random.randint(0, 10) for _ in range(n)]

    def _rand_bool_col(self, n):
        return [random.choice([True, False]) for _ in range(n)]

    def _str_list(self, n, must_have="", count=0):
        s = random.sample(string.ascii_uppercase, k=n)
        if count:
            s = [l + str(count) for l in s]
        if must_have != "":
            s.append(must_have)
        return s


class TestValidationHelpers:

    def _raises_correct_error(self, function, *args,
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

    def dfs_exact_match(self, data, reference):

        if isinstance(data, str):
            data = csvutils.CsvInput(data).read_csv()
        if isinstance(reference, str):
            reference = csvutils.CsvInput(reference).read_csv()

        if set(data.columns) != set(reference.columns):
            return False

        #read csv doesnt precicely read floats
        data = data.round(5)
        reference = reference.round(5)

        return all([reference[col].equals(data[col]) for col in reference.columns])
