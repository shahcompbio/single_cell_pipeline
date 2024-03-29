'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import shutil

import pandas as pd
import pypeliner
from single_cell.utils import csvutils
from single_cell.utils import helpers
from single_cell.utils.singlecell_copynumber_plot_utils import PlotPcolor
from single_cell.workflows.qc_annotation.dtypes import dtypes

from .scripts import classify
from .scripts import fastqscreen_classify
from .scripts import generate_qc


def _get_col_data(df, organism):
    return df['fastqscreen_{}'.format(organism)] - df['fastqscreen_{}_multihit'.format(organism)]


def add_contamination_status(
        infile, outfile, genome_labels,
        reference='grch37', threshold=0.05
):
    data = csvutils.read_csv_and_yaml(infile)

    data = data.set_index('cell_id', drop=False)

    if reference not in genome_labels:
        raise Exception("Could not find the fastq screen counts")

    alts = [col for col in genome_labels if not col == reference]

    data['is_contaminated'] = False

    for altcol in alts:
        perc_alt = _get_col_data(data, altcol) / data['total_reads']
        data.loc[perc_alt > threshold, 'is_contaminated'] = True

    data['is_contaminated'] = data['is_contaminated'].astype('bool')
    csvutils.write_dataframe_to_csv_and_yaml(
        data, outfile, dtypes(genome_labels)
    )


def add_quality(hmmcopy_metrics, alignment_metrics, output, training_data, tempdir, genome_labels):
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

    csvutils.rewrite_csv_file(intermediate_output, output, dtypes=dtypes(genome_labels))


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
    )


def generate_qc_report(
        tempdir, reference_gc, fastqscreen_training_data,
        metrics_df, gc_metrics_df, qc_report, metrics_df_annotated, genome_labels
):
    helpers.makedirs(tempdir)
    fastqscreen_classify.classify_fastqscreen(
        fastqscreen_training_data, metrics_df, metrics_df_annotated, dtypes(genome_labels)
    )

    generate_qc.generate_html_report(tempdir, qc_report, reference_gc, metrics_df_annotated, gc_metrics_df)


def cell_cycle_classifier(hmmcopy_reads, hmmcopy_metrics, alignment_metrics, output, tempdir, genome_labels):
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

    pypeliner.commandline.execute(*cmd)

    cell_cycle_df = pd.read_csv(temp_output)

    cols_cell_cycle = cell_cycle_df.columns.values

    hmm_metrics_df = csvutils.read_csv_and_yaml(hmmcopy_metrics)

    hmm_metrics_df = hmm_metrics_df.merge(cell_cycle_df, on=['cell_id'], how='outer')

    out_dtypes = dtypes(genome_labels)
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
