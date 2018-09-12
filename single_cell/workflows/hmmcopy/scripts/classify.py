'''
Created on Jun 26, 2018

@author: dgrewal
'''
import argparse
import pandas as pd
import numpy as np

import shutil

from sklearn.ensemble import RandomForestClassifier


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
        hmmcopy_data = read_from_h5(hmmcopy_filename, hmmcopy_table)
        alignment_data = read_from_h5(alignment_filename, alignment_table)

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
    predictions = model.predict_proba(data)

    index_good_proba = np.where(model.classes_ == 1)

    assert len(index_good_proba) == 1

    index_good_proba = index_good_proba[0]

    predictions = predictions[:, index_good_proba]

    predictions = dict(zip(data.index, predictions))

    return predictions


def write_to_output(hmmcopy_filename, hmmcopy_tablename, output, predictions):

    data = read_from_h5(hmmcopy_filename, hmmcopy_tablename)

    data['quality'] = data['cell_id'].map(predictions)
    data.quality = data.quality.astype(float)

    with pd.HDFStore(output, 'a', complevel=9, complib='blosc') as outstore:
        outstore.put(hmmcopy_tablename, data, format='table')


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
