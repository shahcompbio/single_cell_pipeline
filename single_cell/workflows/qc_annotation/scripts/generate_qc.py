from __future__ import division
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')

from matplotlib import pyplot as plt
import seaborn as sns

import os

import base64

sns.set(context='talk',
        style='darkgrid',
        # font='Helvetica',
        rc={'axes.titlesize': 12,
            'axes.labelsize': 12,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 12,
            'font.size': 12})


def encode_as_base64(filepath):
    with open(filepath, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read())
    return encoded_string

def load_data(infile, squeeze=False, gc=False):
    df = pd.read_csv(infile, squeeze=squeeze)

    if gc:
        df.index = df.cell_id
        del df['cell_id']
        df.columns = map(int, df.columns.values)

    return df

def flag_outliers(df, upper_limit, lower_limit):
    if df["quality"] < upper_limit and df["quality"] > lower_limit:
        return False
    else:
        return True


def get_row_and_col(df):
    sample, library, row, col = df["cell_id"].split("-")
    df["row"] = row
    df["col"] = col
    return df


def get_spatial_effects(df):
    """
    """
    df = df.apply(get_row_and_col, axis=1)
    df_rows = df.groupby(by="row").mean()
    df_cols = df.groupby(by="col").mean()

    q75_row, q25_row = np.percentile(df_rows, [75, 25])
    q75_col, q25_col = np.percentile(df_cols, [75, 25])
    iqr_row = q75_row - q25_row
    iqr_col = q75_col - q25_col

    lower_limit_rows = q75_row - 1.5 * iqr_row
    upper_limit_rows = q75_row + 1.5 * iqr_row

    lower_limit_cols = q75_col - 1.5 * iqr_col
    upper_limit_cols = q75_col + 1.5 * iqr_col

    df_rows = df_rows.apply(flag_outliers, args=(upper_limit_rows, lower_limit_rows,), axis=1)
    df_cols = df_cols.apply(flag_outliers, args=(upper_limit_cols, lower_limit_cols,), axis=1)

    dropped_rows = df_rows[df_rows]
    dropped_cols = df_cols[df_cols]

    heatmap_df = df.pivot(index="row", columns="col", values="quality")

    return heatmap_df, dropped_rows, dropped_cols


def get_fraction_unmapped(df):
    df["fraction_unmapped"] = df["unmapped_reads"] / df["total_reads"]
    return df


def get_quality_pass_metrics(df):
    cells_pass_df = df[df["quality"] >= 0.75]
    # Calculate total cells pass

    # Calculate cells pass by experimental condition and by cell call
    cells_pass = cells_pass_df.groupby(["experimental_condition", "cell_call"]).count()["cell_id"]
    # Calculate median coverage depth
    median_coverage_depth = cells_pass_df.groupby(["experimental_condition", "cell_call"]).median()["coverage_depth"]

    return cells_pass, median_coverage_depth


def get_dropout_metrics(df):
    cells_dropout_df = df[df["total_mapped_reads"] < 250000]
    # Calculate cells dropout
    cells_dropout = cells_dropout_df.groupby(["experimental_condition", "cell_call"]).count()["cell_id"]

    # Calculate total breakdown
    total_breakdown = df.groupby(["experimental_condition", "cell_call"]).count()["cell_id"]

    # Calculate fraction cells dropout
    fraction_cells_dropout = cells_dropout.divide(total_breakdown, fill_value=0)

    return cells_dropout, fraction_cells_dropout


def get_non_dropout_metrics(df):
    cells_non_dropout_df = df[df["total_mapped_reads"] >= 250000]

    # Calculate cells flagged
    all_flagged_cells = cells_non_dropout_df[cells_non_dropout_df["quality"] < 0.75]

    cells_flagged = all_flagged_cells.groupby(["experimental_condition", "cell_call"]).count()["cell_id"]

    # Calculate median quality of non-dropouts
    median_quality = cells_non_dropout_df.groupby(["experimental_condition", "cell_call"]).median()["quality"]

    # Calculate total breakdown
    total_breakdown = df.groupby(["experimental_condition", "cell_call"]).count()["cell_id"]

    # Calculate fraction of cells flagged
    fraction_cells_flagged = cells_flagged.divide(total_breakdown, fill_value=0)

    return cells_flagged, median_quality, fraction_cells_flagged


def get_cells_ref_species(df):
    cells_not_other_species_df = df[df["total_reads"] < 250000]

    if cells_not_other_species_df.empty:
        return pd.Series()

    # Calculate percent contamination
    cells_not_other_species = cells_not_other_species_df.apply(get_fraction_unmapped, axis=1)
    cells_not_other_species = cells_not_other_species[cells_not_other_species["fraction_unmapped"] < 0.5]
    median_fraction_background = cells_not_other_species.groupby(["experimental_condition", "cell_call"]).median()[
        "fraction_unmapped"]

    return median_fraction_background


def get_cells_other_species_df(df):
    # Calculate total breakdown
    total_breakdown = df.groupby(["experimental_condition", "cell_call"]).count()["cell_id"]

    cells_other_species_df = df[df["total_reads"] >= 250000]
    unmapped_df = cells_other_species_df.apply(get_fraction_unmapped, axis=1)
    cells_other_species = \
        unmapped_df[unmapped_df["fraction_unmapped"] >= 0.5].groupby(["experimental_condition", "cell_call"]).count()[
            "cell_id"]

    # Calculate percent cells other species
    fraction_cells_other_species = cells_other_species.divide(total_breakdown, fill_value=0)

    return cells_other_species, fraction_cells_other_species


def generate_contamination_qc_table(df):
    median_fraction_bg = get_cells_ref_species(df)
    cells_other_species, fraction_cells_other_species = get_cells_other_species_df(df)

    metrics = [
        cells_other_species.rename("cells_other_species"),
        fraction_cells_other_species.rename("fraction_cells_other_species"),
        median_fraction_bg.rename("median_fraction_background")
    ]

    # Create the output table
    metrics = pd.concat(metrics, axis=1)

    return metrics


def generate_quality_qc_table(df):
    cells_pass, median_cov_depth = get_quality_pass_metrics(df)
    cells_dropout, fraction_cells_dropout = get_dropout_metrics(df)
    cells_flagged, median_quality, fraction_cells_flagged = get_non_dropout_metrics(df)

    # Calculate total breakdown
    total_breakdown = df.groupby(["experimental_condition", "cell_call"]).count()["cell_id"]
    # Calculate percentage of cells pass
    fraction_cell_pass = cells_pass.divide(total_breakdown, fill_value=0)
    # Calculate median reads
    median_reads = df.groupby(["experimental_condition", "cell_call"]).median()["total_reads"]

    metrics = [
        cells_pass.rename("cells_pass"),
        total_breakdown.rename("total_breakdown"),
        fraction_cell_pass.rename("fraction_cell_pass"),
        cells_dropout.rename("cells_dropout"),
        fraction_cells_dropout.rename("fraction_cells_dropout"),
        cells_flagged.rename("cells_flagged"),
        fraction_cells_flagged.rename("fraction_cells_flagged"),
        median_quality.rename("median_quality"),
        median_reads.rename("median_reads"),
        median_cov_depth.rename("median_coverage_depth"),
    ]
    # Create the output table
    metrics = pd.concat(metrics, axis=1)

    return metrics


def generate_library_metrics(df, gc_data, reference_gc):
    cells_pass_df = df[df["quality"] >= 0.75]
    # Calculate total cells pass
    total_cells_pass = cells_pass_df.count()["cell_id"]

    # Find the error
    gc_curve = get_gc_curve(gc_data)
    error_curve = reference_gc.sub(gc_curve).abs()
    sum_of_diffs = float(error_curve.sum())

    df = pd.DataFrame.from_dict(
        {
            'Cell Yield': int(total_cells_pass),
            'sum of GC bias differences': sum_of_diffs,
        },
        orient='index'
    )

    return df


def get_gc_curve(df):
    """
    """
    df_average = df.mean().to_frame(name="average")

    q75, q25 = np.percentile(df_average["average"], [75, 25])
    iqr = q75 - q25

    # Remove extreme outliers
    upper_limit = q75 + 3 * iqr
    lower_limit = q25 - 3 * iqr
    average = df_average.iloc[0:75].average.apply(remove_outliers, args=(upper_limit, lower_limit,))
    average = average.fillna(0)

    return average


def remove_outliers(average, upper_limit, lower_limit):
    if average < upper_limit and average > lower_limit:
        return average
    else:
        return np.nan


def plot_gc_curve(gc_data, reference_gc, output):
    """
    """

    gc_curve = get_gc_curve(gc_data)

    reference_gc.index = list(reference_gc.index)
    gc_curve.index = list(gc_curve.index)

    # Compare with reference gc curve
    reference_gc.index = list(reference_gc.index)
    gc_curve.index = list(gc_curve.index.astype(int))

    plt.grid(b=True, which='both', axis='both')
    plt.plot(reference_gc, color='b', label='reference curve')
    plt.plot(gc_curve, 'o', mfc='none', mec='r', mew=1, label='Normalized Coverage')
    plt.title("GC Bias Plot")
    plt.xlabel("GC% of 100 base windows")
    plt.ylabel("Fraction of normalized coverage")
    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig(output)
    plt.clf()


def plot_heatmap(data, output):
    """
    """
    heatmap_df, dropped_rows, dropped_cols = get_spatial_effects(data[["cell_id", "quality"]])

    mask = heatmap_df.isna()
    # sns.set(font_scale=0.5)
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(
        heatmap_df,
        mask=mask,
        xticklabels=True,
        yticklabels=True,
        cbar=True,
        cbar_kws={"shrink": .3}
    )
    for tick in ax.get_xticklabels():
        if tick.get_text() in dropped_cols:
            tick.set_color("r")

    for tick in ax.get_yticklabels():
        if tick.get_text() in dropped_rows:
            tick.set_color("r")
    plt.xlabel("Column")
    plt.ylabel("Row")
    plt.title("Quality per Cell", fontsize=12)
    plt.tight_layout()
    plt.savefig(output)
    plt.clf()


def pretty():
    def hover(hover_color="#e8e8e8"):
        return dict(selector="tr:hover",
                    props=[("background-color", "%s" % hover_color)])

    return [
        hover(),
        dict(selector="th", props=[('font-size', '11pt'), ("font-family", "Helvetica"),
                                   ("color", 'black'),
                                   ('background-color', 'rgb(232, 232, 232)'),
                                   ("text-align", "right")]),
        dict(selector="td", props=[("font-family", "Helvetica"),
                                   ("text-align", "right")]),
        dict(selector="tr", props=[("line-height", "11px")]),
        dict(selector="caption", props=[("caption-side", "bottom")])
    ]


def generate_html(dataframes, pngs, html_file):
    html_elements = []

    for header, df in dataframes:
        html_elements.append("<h3>{}</h3>\n".format(header))
        html_elements.append(df.style.set_table_styles(pretty()).render())

    for header, image in pngs:
        image = encode_as_base64(image)
        html_elements.append("<h3>{}</h3>\n".format(header))
        html_elements.append('<img src="data:image/png;base64,{}" alt="beastie.png" scale="0">'.format(image))

    with open(html_file, 'w') as html_out:
        for element in html_elements:
            html_out.write(element)
            html_out.write('\n')


def generate_html_report(tempdir, html, reference_gc, metrics, gc_metrics):
    gc_plot = os.path.join(tempdir, "gc.png")
    heatmap = os.path.join(tempdir, 'heatmap.png')

    reference_gc = load_data(reference_gc, squeeze=True)

    data = load_data(metrics)
    gc_data = load_data(gc_metrics, gc=True)

    quality_df = generate_quality_qc_table(data)
    contamination_df = generate_contamination_qc_table(data)
    library_df = generate_library_metrics(data, gc_data, reference_gc)

    plot_gc_curve(gc_data, reference_gc, gc_plot)
    plot_heatmap(data, heatmap)

    generate_html(
        [
            ('Metrics', quality_df),
            ('Contamination', contamination_df),
            ('Other Metrics', library_df),
        ],
        [
            ('Heatmap', heatmap),
            ('GC plot', gc_plot)
        ],
        html
    )