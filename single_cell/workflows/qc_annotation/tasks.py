'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner

from single_cell.utils import csvutils
from single_cell.utils import helpers

from scripts import classify
from scripts import generate_qc


def annotate_metrics(
        metrics, output, sample_info, cells):
    """
    adds sample information to metrics in place
    """

    metrics = csvutils.read_csv_and_yaml(metrics)

    for cellid in cells:
        cellinfo = sample_info[cellid]

        for colname, value in cellinfo.iteritems():
            metrics.loc[metrics["cell_id"] == cellid, colname] = value

    csvutils.write_dataframe_to_csv_and_yaml(metrics, output)


def add_quality(hmmcopy_metrics, alignment_metrics, output, training_data, tempdir):
    helpers.makedirs(tempdir)

    intermediate_output = os.path.join(tempdir, 'metrics_with_quality.csv')

    model = classify.train_classifier(training_data)

    feature_names = model.feature_names_

    data = classify.load_data(hmmcopy_metrics, alignment_metrics,
                              feature_names)

    predictions = classify.classify(model, data)

    classify.write_to_output(
        hmmcopy_metrics,
        intermediate_output,
        predictions)

    csvutils.prep_csv_files(intermediate_output, output)


def merge_metrics(hmmcopy_metrics, alignment_metrics, merged_output):
    csvutils.merge_csv(
        [
            hmmcopy_metrics,
            alignment_metrics
        ],
        merged_output,
        'outer',
        ['cell_id'],
        header=True
    )


def generate_qc_report(tempdir, reference_gc, metrics_df, gc_metrics_df, qc_report):
    helpers.makedirs(tempdir)

    generate_qc.generate_html_report(tempdir, qc_report, reference_gc, metrics_df, gc_metrics_df)


def cell_cycle_classifier(hmmcopy_reads, hmmcopy_metrics, alignment_metrics, output, docker_image=None):
    cmd = [
        'cell_cycle_classifier',
        'train-classify',
        hmmcopy_reads,
        hmmcopy_metrics,
        alignment_metrics,
        output
    ]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)
