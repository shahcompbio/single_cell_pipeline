'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import shutil

import pandas as pd
import pypeliner
from ete3 import Tree
from single_cell.utils import csvutils
from single_cell.utils import helpers
from single_cell.utils.singlecell_copynumber_plot_utils import PlotPcolor
from single_cell.workflows.qc_annotation.dtypes import dtypes

from .scripts import classify
from .scripts import generate_qc


def add_corrupt_tree_order(corrupt_tree, metrics, output):
    """
    adds corrupt tree order to metrics
    """

    with open(corrupt_tree) as newickfile:
        newickdata = newickfile.readline()
        assert newickfile.readline() == ''

    tree = Tree(newickdata, format=1)

    leaves = [node.name for node in tree.traverse("levelorder")]
    leaves = [val[len('cell_'):] for val in leaves if val.startswith("cell_")]

    ordering = {val: i for i, val in enumerate(leaves)}

    metrics = csvutils.read_csv_and_yaml(metrics)

    cells = metrics.cell_id

    for cellid in cells:
        order = ordering.get(cellid, float('nan'))
        metrics.loc[metrics["cell_id"] == cellid, "order_corrupt_tree"] = order

    col_dtype = dtypes()['metrics']['order_corrupt_tree']
    metrics['order_corrupt_tree'] = metrics['order_corrupt_tree'].astype(col_dtype)

    csvutils.write_dataframe_to_csv_and_yaml(
        metrics, output, dtypes()['metrics'], write_header=True
    )


def annotate_metrics(
        metrics, output, sample_info, cells):
    """
    adds sample information to metrics in place
    """

    metrics = csvutils.read_csv_and_yaml(metrics)

    for cellid in cells:
        cellinfo = sample_info[cellid]

        for colname, value in cellinfo.items():
            metrics.loc[metrics["cell_id"] == cellid, colname] = value

    csvutils.write_dataframe_to_csv_and_yaml(metrics, output, dtypes()['metrics'])


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

    csvutils.prep_csv_files(intermediate_output, output, dtypes=dtypes()['metrics'])


def merge_metrics(hmmcopy_metrics, alignment_metrics, merged_output):
    csvutils.merge_csv(
        [
            hmmcopy_metrics,
            alignment_metrics
        ],
        merged_output,
        'outer',
        ['cell_id'],
        write_header=False,
        dtypes=dtypes()['metrics']
    )


def generate_qc_report(tempdir, reference_gc, metrics_df, gc_metrics_df, qc_report):
    helpers.makedirs(tempdir)

    generate_qc.generate_html_report(tempdir, qc_report, reference_gc, metrics_df, gc_metrics_df)


def cell_cycle_classifier(hmmcopy_reads, hmmcopy_metrics, alignment_metrics, output, tempdir, docker_image=None):
    helpers.makedirs(tempdir)
    temp_output = os.path.join(tempdir, 'cell_cycle_output.csv')

    cmd = [
        'cell_cycle_classifier',
        'train-classify',
        hmmcopy_reads,
        hmmcopy_metrics,
        alignment_metrics,
        temp_output
    ]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    cell_cycle_df = pd.read_csv(temp_output)

    cols_cell_cycle = cell_cycle_df.columns.values

    hmm_metrics_df = csvutils.read_csv_and_yaml(hmmcopy_metrics)

    hmm_metrics_df = hmm_metrics_df.merge(cell_cycle_df, on=['cell_id'], how='outer')

    out_dtypes = dtypes()['metrics']
    for colname in cols_cell_cycle:
        hmm_metrics_df[colname] = hmm_metrics_df[colname].astype(out_dtypes[colname])

    csvutils.write_dataframe_to_csv_and_yaml(hmm_metrics_df, output, out_dtypes)


def filter_plot_tar(metrics, src_tar, pass_tar, fail_tar, tempdir, filters):
    allplots = os.path.join(tempdir, 'allplots')
    helpers.makedirs(allplots)
    helpers.extract_tar(src_tar, allplots)

    metrics_data = csvutils.read_csv_and_yaml(metrics)
    all_cells = metrics_data.cell_id.tolist()

    metrics_data = helpers.filter_metrics(metrics_data, filters)
    good_cells = metrics_data.cell_id.tolist()
    bad_cells = [cell for cell in all_cells if cell not in good_cells]

    plotdir = os.path.join(tempdir, 'segs_pass')
    helpers.makedirs(plotdir)
    for cell in good_cells:
        src_path = os.path.join(allplots, 'segments', '{}_{}.png'.format(cell, 'segments'))
        dest_path = os.path.join(plotdir, '{}_{}.png'.format(cell, 'segments'))
        shutil.copyfile(src_path, dest_path)
    helpers.make_tarfile(pass_tar, plotdir)

    plotdir = os.path.join(tempdir, 'segs_fail')
    helpers.makedirs(plotdir)
    for cell in bad_cells:
        src_path = os.path.join(allplots, 'segments', '{}_{}.png'.format(cell, 'segments'))
        dest_path = os.path.join(plotdir, '{}_{}.png'.format(cell, 'segments'))
        shutil.copyfile(src_path, dest_path)
    helpers.make_tarfile(fail_tar, plotdir)


def get_good_cells(metrics, cell_filters):
    metrics_data = csvutils.read_csv_and_yaml(metrics)

    if not cell_filters:
        return metrics_data.cell_id.tolist()

    metrics_data = helpers.filter_metrics(metrics_data, cell_filters)

    return metrics_data.cell_id.tolist()


def plot_pcolor(infile, metrics, output, corrupt_tree=None, plot_title=None,
                column_name=None, plot_by_col=None,
                chromosomes=None, max_cn=None,
                scale_by_cells=None, color_by_col=None,
                mappability_threshold=None, cell_filters=None
                ):
    cells = get_good_cells(metrics, cell_filters)

    plot = PlotPcolor(
        infile, metrics, output, plot_title=plot_title,
        column_name=column_name, plot_by_col=plot_by_col,
        chromosomes=chromosomes,
        max_cn=max_cn,
        scale_by_cells=scale_by_cells,
        color_by_col=color_by_col,
        corrupt_tree=corrupt_tree,
        mappability_threshold=mappability_threshold,
        cells=cells,
    )
    plot.main()
