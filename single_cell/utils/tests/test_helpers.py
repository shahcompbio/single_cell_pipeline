import string
import single_cell.utils.csvutils as csvutils
import pandas as pd
import random
import os
import yaml as y


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


class AnnotationHelpers(TestInputs, TestValidationHelpers):
    """
    helper functions for testing the annotate_csv csvutils function
    """

    def sim_ann_input(self, annotate_col, dtypes):
        """
        make annotation dict for test annotate_csv.
        :param annotate_col: column to create annotation on
        :param dtypes: dtypes for annotation dict
        """
        annotations = {cell_id: {col: self.simulate_col(dtype, 1)[0]
                                 for col, dtype in dtypes.items()}
                       for cell_id in annotate_col}

        return annotations

    def make_ann_test_inputs(self, temp, length, dtypes, ann_dtypes,
                             write_header=True, on="cell_id"):
        """
        make inputs to test annotate_csv
        :param temp: tempdir to test in
        :param length: length of test dfs
        :param dtypes: dtypes of test dfs
        :param ann_dtypes: dtypes of test annotation dict
        :param write_header: T/F write header post-annotation
        :param on: col to annotate on
        """
        df = self.make_test_dfs([dtypes], length)
        csv = self.write_dfs(temp, df, [dtypes], write_header)

        df = df[0]
        csv = csv[0]

        annotation_input = self.sim_ann_input(df[on].tolist(), ann_dtypes)

        return csv, annotation_input

    def base_annotation_test(self, temp, length, dtypes, ann_dtypes, head=True,
                             on="cell_id"):
        """
        base test of annotate_csv
        :param temp: tempdir to test in
        :param length: length of test dfs
        :param dtypes: dtypes of test dfs
        :param ann_dtypes: dtypes of test annotation dict
        :param write_header: T/F write header post-annotation
        :param on: col to annotate on:
        """
        csv, annotation = self.make_ann_test_inputs(temp, length,
                                                    dtypes, ann_dtypes,
                                                    write_header=head,
                                                    on=on)

        annotated = os.path.join(temp, "annotated.csv.gz")

        csvutils.annotate_csv(csv, annotation, annotated, ann_dtypes,
                              write_header=head, on=on)

        return csv, annotation, annotated

    def validate_annotation_test(self, csv, annotation, annotated, on):
        """
        test that annotate_csv produced correct result in tests
        :param csv: csv that was annotated
        :param annotation: input annotation dict
        :param annotated: csv post-annotation
        :param on: col on which csv was annotated
        :return:
        """
        if isinstance(csv, str):
            csv = csvutils.CsvInput(csv).read_csv()

        annotation_comparable = {on: list(annotation.keys())}
        annotation_as_list = list(annotation.values())

        for k in annotation_as_list[0].keys():
            annotation_comparable[k] = [ann[k] for ann in annotation_as_list]

        annotation_comparable = pd.DataFrame(annotation_comparable)

        compare = csv.merge(annotation_comparable, how="outer", on=on)

        return os.path.exists(annotated) \
               and self.dfs_exact_match(compare, annotated)



class ConcatHelpers(TestInputs, TestValidationHelpers):
    """
    helpers for testing concat functions
    """
    def base_test_concat(self, length, dtypes, write=False, dir=None,
                         get_ref=False, write_head=True):
        """
        base concatenation test; make dfs, write them,
        get a "reference" concat output
        :param length: number of rows in df; use length not n_rows so as to not use fixture
        :param dtypes: list of dtypes to use in creation of test dfs
        :param write: T/F :write out dfs to csvs
        :param dir: directory to write to if writing
        :param get_ref: T/F: get excpected output from concatenation using
        pandas built-ins
        :param write_head: T/F: write header when writing out to csvs
        :return: dfs and or csvs and or expected output
        """

        dfs = self.make_test_dfs(dtypes, length)

        outs = [dfs]

        if write:
            csvs = self.write_dfs(dir, dfs, dtypes, write_head)
            outs.append(csvs)

        if get_ref:
            ref = pd.concat(dfs, ignore_index=True)
            outs.append(ref)

        return tuple(outs)


class WriteHelpers(TestInputs, TestValidationHelpers):
    """
    helpers class for testing of csvutils writing
    """
    def write_file_successful(self, df, csv, wrote_header=True):
        head = "infer"
        if not wrote_header:
            head = None
        csv = pd.read_csv(csv, sep=",", header=head)

        if not wrote_header:
            csv.columns = df.columns

        if not df.equals(csv):
            return False
        return True

    def metadata_write_successful(self, dtypes, yaml_file):
        """
        validate writing of dtypes to yaml
        :param dtypes: written dtypes
        :param yaml_file: output .yaml
        :return:
        """
        yaml_loader = y.load(open(yaml_file), Loader=y.FullLoader)["columns"]
        yaml = {}

        for dtype in yaml_loader:
            yaml[dtype["name"]] = dtype["dtype"]

        if yaml != dtypes:
            return False
        return True


class MergeHelpers(TestInputs, TestValidationHelpers):
    """
    helpers class for testing of merging
    """
    def make_mergeable_test_dfs(self, dtypes, shared, length):
        """
        make dfs that can be merged
        :param dtypes: list of dtypes corresponding to test dfs
        :param shared: columns shared by dfs, should all be in dtypes
        :param length: length of test dfs
        :return:
        """
        for dtype_set in dtypes:
            assert set(shared).issubset(set(dtype_set.keys()))

        dfs = self.make_test_dfs(dtypes, length)

        shared_values = dfs[0][shared]

        for df in dfs:
            df.update(shared_values)

        return dfs

    def base_merge_test(self, length, how, on, suffixes, dtypes, write=False,
                        get_ref=True, dir=None, write_head=True):
        """
        base test for all merge fnctions
        :param length: length of test dfs
        :param how: how to merge (i.e. direction)
        :param on: col(s) to merge on
        :param suffixes: suffixes for merged df (see pandas docs)
        :param dtypes: list of dtypes for creation of test dfs
        :param write: T/F write out test dfs before merging
        :param get_ref: get expected merge output
        :param dir: temp dir to test in
        :param write_head: T/F write header to csv (before and after testing)
        """
        dfs = self.make_mergeable_test_dfs(dtypes, on, length)

        outs = [dfs]
        if not write and not get_ref:
            return dfs

        if write:
            csvs = self.write_dfs(dir, dfs, dtypes, write_head)
            outs.append(csvs)

        if get_ref:
            ref = dfs[0].merge(dfs[1], how=how, on=on, suffixes=suffixes)
            outs.append(ref)

        return tuple(outs)

