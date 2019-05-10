'''
Created on Jun 26, 2018

@author: dgrewal
'''
import argparse
import pandas as pd
import numpy as np
import single_cell.utils.helpers
import shutil

from sklearn.ensemble import RandomForestClassifier

from single_cell.utils import csvutils


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--training_data',
                        help='specify the path to training data h5')

    parser.add_argument('--alignment_metrics',
                        help='specify path to the alignment metrics h5')

    parser.add_argument('--hmmcopy_metrics',
                        help='specify path to the hmmcopy metrics h5')

    parser.add_argument('--output',
                        help='specify output file'
                        )

    parser.add_argument('--alignment_table',
                        help='table names in h5 for alignment'
                        )

    parser.add_argument('--hmmcopy_tables',
                        nargs='*',
                        help='table names in h5 for hmmcopy'
                        )

    args = parser.parse_args()

    return args


def read_from_h5(filename, tablename):
    with pd.HDFStore(filename) as h5store:
        data = h5store[tablename]
    return data


# def read_from_csv(filename, gzipped=False):
#     compression = 'gzip' if gzipped else None
#     data = pd.read_csv(filename, compression=compression)
#     return data


def read_data(filename, tablename, gzipped=True):
    fileformat = single_cell.utils.helpers.get_file_format(filename)

    if fileformat == 'h5':
        data = read_from_h5(filename, tablename)
    elif fileformat == 'csv':
        data = csvutils.read_csv_and_yaml(filename)
    elif fileformat == 'gzip':
        data = csvutils.read_csv_and_yaml(filename)
    else:
        raise Exception("unknown file format")

    return data


def train_classifier(filename):
    training_data = read_from_h5(filename, '/training_data')

    labels = training_data["label"]

    del training_data["label"]

    clf = RandomForestClassifier(n_estimators=500)

    model = clf.fit(training_data, labels)

    features = training_data.columns.values
    model.feature_names_ = features

    return model


def load_data(hmmcopy_filename, alignment_filename,
              hmmcopy_tables, alignment_table, colnames):
    for hmmcopy_table in hmmcopy_tables:
        hmmcopy_data = read_data(hmmcopy_filename, hmmcopy_table)
        alignment_data = read_data(alignment_filename, alignment_table)

        hmmcopy_data = hmmcopy_data.set_index('cell_id')
        alignment_data = alignment_data.set_index('cell_id')

        data = []
        for colname in colnames:
            if colname in hmmcopy_data:
                coldata = hmmcopy_data[colname]
            else:
                coldata = alignment_data[colname]

            if colname == 'scaled_halfiness':
                # haploid poison adds inf, replace with big number since 0 is considered good
                # and we want to score to decrease
                coldata = coldata.replace(np.inf, 1e10)
            data.append(coldata)

        data = pd.concat(data, axis=1)

        data = data.replace(-np.inf, np.nan)
        data = data.fillna(0)

        yield (hmmcopy_table, data)


def classify(model, data):
    ##TODO: remove this with v0.2.3, also handled in collect_metrics
    # picardtools sometimes reports missing as ?
    data = data.replace('?', 0)
    predictions = model.predict_proba(data)

    index_good_proba = np.where(model.classes_ == 1)

    assert len(index_good_proba) == 1

    index_good_proba = index_good_proba[0]

    predictions = predictions[:, index_good_proba]

    predictions = dict(zip(data.index, predictions))

    return predictions


def write_to_hdf(output, hmmcopy_tablename, data):
    with pd.HDFStore(output, 'a', complevel=9, complib='blosc') as outstore:
        outstore.put(hmmcopy_tablename, data, format='table')


def write_to_csv(output, data, gzipped=False):
    compression = 'gzip' if gzipped else None
    data.to_csv(output, index=False, compression=compression)


def write_to_output(hmmcopy_filename, hmmcopy_tablename, output, predictions):
    data = read_data(hmmcopy_filename, hmmcopy_tablename)

    data['quality'] = data['cell_id'].map(predictions)
    data.quality = data.quality.astype(float)

    fileformat = single_cell.utils.helpers.get_file_format(output)

    if fileformat == 'h5':
        write_to_hdf(output, hmmcopy_tablename, data)
    elif fileformat == 'csv':
        write_to_csv(output, data)
    elif fileformat == "gzip":
        write_to_csv(output, data, gzipped=True)
    else:
        raise Exception("unknown file format")


if __name__ == "__main__":
    args = parse_args()

    shutil.copy(args.hmmcopy_metrics, args.output)

    model = train_classifier(args.training_data)

    feature_names = model.feature_names_

    data = load_data(args.hmmcopy_metrics, args.alignment_metrics,
                     args.hmmcopy_tables, args.alignment_table,
                     feature_names)

    for hmmcopy_table, tabledata in data:
        predictions = classify(model, tabledata)

        write_to_output(
            args.hmmcopy_metrics,
            hmmcopy_table,
            args.output,
            predictions)
