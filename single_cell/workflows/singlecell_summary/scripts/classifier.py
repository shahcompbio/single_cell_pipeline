'''
Created on Mar 22, 2018

@author: tchan
'''

import argparse
import pandas as pd
from sklearn.externals import joblib


def read_metrics(path):

    df = pd.read_csv(path,
            na_values=['', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN',
                       '-NaN', '-nan', '1.#IND', '1.#QNAN', 'N/A', 'NA', 'NULL',
                        'NaN', 'n/a', 'nan', 'null', '?'])

    return df


def load_model(path):
    return joblib.load(path)


def cleanup_metrics(metrics):
    drop_columns = ['all_heatmap_order', 'sample_well', 'i5_barcode', 'i7_barcode', 'sample_type',
                    'sample_plate', 'cell_call', 'cell_id', 'experimental_condition']

    metrics = metrics.drop(columns=drop_columns)

    metrics = metrics.fillna(0)

    return metrics


def classify(data, model):

    model_features = model.feature_names

    rm_cols = [col for col in data.columns.values if col not in model_features]

    if rm_cols:
        data = data.drop(columns=rm_cols)

    data = data[model_features]

    assert list(data.columns.values) == model.feature_names

    predictions = model.predict_proba(data)

    predictions = [v[1] for v in predictions]

    return predictions


def predict(metrics, model):

    cells = metrics.cell_id

    metrics = cleanup_metrics(metrics)

    predictions = classify(metrics, model)

    df = pd.DataFrame()
    df["cell_id"] = cells
    df["probability_good"] = predictions
    return df


def add_column_to_metrics(predictions, metrics):

    metrics = metrics.merge(predictions, how="inner", on="cell_id")

    return metrics


def write(df, output):

    df.to_csv(output, index=False, na_rep="NA")


def main(metrics, model, output):

    metrics = read_metrics(metrics)

    model = load_model(model)

    predictions = predict(metrics, model)

    metrics = add_column_to_metrics(predictions, metrics)

    write(metrics, output)


def get_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("metrics_csv")

    parser.add_argument("model")

    parser.add_argument("destination")

    args = vars(parser.parse_args())
    return args


if __name__ == '__main__':
    args = get_args()
    main(args['metrics_csv'], args['model'], args['destination'])
