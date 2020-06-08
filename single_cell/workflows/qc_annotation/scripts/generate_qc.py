from __future__ import division

import matplotlib
import numpy as np
import pandas as pd
matplotlib.use('Agg')

from matplotlib import pyplot as plt
import seaborn as sns

import os

import base64

from single_cell.utils import csvutils

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
    return encoded_string.decode("utf-8")


def load_reference(infile):
    df = pd.read_csv(infile, squeeze=True)
    return df


def load_data(infile, gc=False):
    df = csvutils.read_csv_and_yaml(infile)

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
    information = df["cell_id"].split("-")
    if len(information) == 4:
        sample, library, row, col = information
    else:
        sample, library, row, col = information[0], "_".join(information[1:-2]) ,information[-2], information[-1]
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
    df["fraction_unmapped"] = np.divide(df["unmapped_reads"], df["total_reads"])
    return df

def get_dropout_metrics(dropout_df, total_breakdown):
    # Calculate cells dropout
    cells_dropout = dropout_df.groupby(["Experimental Condition", "Cell Call"]).count()["cell_id"]

    # Calculate fraction cells dropout
    fraction_cells_dropout = cells_dropout.divide(total_breakdown, fill_value=0)

    metrics = {
        "cells_dropout": cells_dropout,
        "fraction_cells_dropout": fraction_cells_dropout
    }

    return metrics

def get_non_dropout_metrics(non_dropout_df, total_breakdown):
    # Calculate cells flagged (Cells with >250k reads, but <0.75 quality (LQ cells))
    all_flagged_cells = non_dropout_df[non_dropout_df["quality"] < 0.75]
    count_cells_flagged = all_flagged_cells.groupby(["Experimental Condition", "Cell Call"]).count()["cell_id"]

    # Calculate fraction of cells flagged
    fraction_cells_flagged = count_cells_flagged.divide(total_breakdown, fill_value=0)

    non_dropout_cells = non_dropout_df.groupby(["Experimental Condition", "Cell Call"]).count()["cell_id"]
    
    human_count = salmon_count = mouse_count = None
    human_ratio = salmon_ratio = mouse_ratio = None
    
    if "species" in non_dropout_df:   
        human_groupped = non_dropout_df[non_dropout_df["species"]=="grch37"].groupby(["Experimental Condition", "Cell Call"])
        mouse_groupped = non_dropout_df[non_dropout_df["species"]=="mm10"].groupby(["Experimental Condition", "Cell Call"])
        salmon_groupped = non_dropout_df[non_dropout_df["species"]=="salmon"].groupby(["Experimental Condition", "Cell Call"])

        human_count = human_groupped.count()["cell_id"]
        mouse_count = mouse_groupped.count()["cell_id"]
        salmon_count = salmon_groupped.count()["cell_id"]

        human_ratio = human_count.divide(non_dropout_cells, fill_value=0)
        mouse_ratio = mouse_count.divide(non_dropout_cells, fill_value=0)
        salmon_ratio = salmon_count.divide(non_dropout_cells, fill_value=0)

    metrics = {
            "count_cells_flagged": count_cells_flagged,
            "fraction_cells_flagged": fraction_cells_flagged,
            "non_dropout_cells": non_dropout_cells,
            "human_count": human_count,
            "mouse_count": mouse_count,
            "salmon_count": salmon_count,
            "human_ratio": human_ratio,
            "mouse_ratio": mouse_ratio,
            "salmon_ratio": salmon_ratio
        }


    return metrics


def get_hq_metrics(hq_df, total_breakdown):
    #get metrics of high quality cells
    #group by experimental condition and cell call
    hq_cell_groupped = hq_df.groupby(["Experimental Condition", "Cell Call"])
    #number of HQ cells
    hq_cell_count = hq_cell_groupped.count()["cell_id"]
    #percentage of HQ cells
    hq_cell_percentage = hq_cell_count.divide(total_breakdown, fill_value=0)
    #mean quality
    hq_mean_quality = hq_cell_groupped.mean()["quality"]
    #median quality
    hq_median_quality = hq_cell_groupped.median()["quality"]
    #median reads
    hq_median_reads = hq_cell_groupped.median()["total_reads"]
    #median coverage depth
    hq_median_coverage_depth = hq_cell_groupped.median()["coverage_depth"]
    #Median % of unmapped reads in HQ cells
    median_unmapped_ratio = hq_cell_groupped.median()["fraction_unmapped"]
    #Median % of unknown reads in HQ cells
    median_unknown_reads_ratio = None
    salmon = None
    human = None
    mouse = None
    if "fastqscreen_nohit_ratio" in hq_df:
        median_unknown_reads_ratio = hq_cell_groupped.median()["fastqscreen_nohit_ratio"]

    if "species" in hq_df:
        #species counts
        salmon = hq_df[hq_df["species"]=="salmon"].groupby(["Experimental Condition", "Cell Call"]).count()["cell_id"]
        human = hq_df[hq_df["species"]=="grch37"].groupby(["Experimental Condition", "Cell Call"]).count()["cell_id"]
        mouse = hq_df[hq_df["species"]=="mm10"].groupby(["Experimental Condition", "Cell Call"]).count()["cell_id"]
    
    hq_metrics = {
            "hq_cell_count": hq_cell_count,
            "hq_cell_percentage": hq_cell_percentage,
            "hq_mean_quality": hq_mean_quality,
            "hq_median_quality": hq_median_quality,
            "hq_median_reads": hq_median_reads,
            "hq_median_coverage_depth": hq_median_coverage_depth,
            "hq_salmon_count": salmon,
            "hq_human_count": human,
            "hq_mouse_count": mouse,
            "median_unmapped_ratio": median_unmapped_ratio,
            "median_unknown_reads_ratio": median_unknown_reads_ratio,
        }

    return hq_metrics

def generate_qc_table(df):
    #assign species labels to the cells
    df = get_fraction_unmapped(df)
    df = df.rename(columns={"experimental_condition": "Experimental Condition", "cell_call": "Cell Call"})
    #get high quality cells
    hq_df = df[df["quality"] >= 0.75]

    #get low quality cells
    lq_df = df[df["quality"] < 0.75]

    #get dropout cells
    dropout_df = df[df["total_mapped_reads"] < 250000]

    #get non-dropout cells
    non_dropout_df = df[df["total_mapped_reads"] >= 250000]

    #total number of cells
    total_breakdown = df.groupby(["Experimental Condition", "Cell Call"]).count()["cell_id"]

    #manual place holder
    place_holder = total_breakdown.copy().apply(lambda x: "place_holder")
    
    ### metrics for high quality cells
    hq_metrics = get_hq_metrics(hq_df, total_breakdown)
    
    ### metrics for dropout
    dropout_metrics = get_dropout_metrics(dropout_df, total_breakdown)

    ### metrics for non-dropout
    non_dropout_metrics = get_non_dropout_metrics(non_dropout_df, total_breakdown)

    metrics_lst = [
        total_breakdown.rename("All cells"),
        place_holder.rename("line_break_0"),

        hq_metrics["hq_cell_count"].rename("HQ cells"),
        hq_metrics["hq_cell_percentage"].astype(float).map("{:.2%}".format).rename("% of HQ cells"),
        place_holder.rename("line_break_1"),

        dropout_metrics["cells_dropout"].rename("Cells with <250k reads (Dropout cells)"),
        dropout_metrics["fraction_cells_dropout"].astype(float).map("{:.2%}".format).rename("% of Dropout cells"),
        place_holder.rename("line_break_2"),

        non_dropout_metrics["count_cells_flagged"].rename("Cells with >250k reads, but <0.75 quality"),
        non_dropout_metrics["fraction_cells_flagged"].astype(float).map("{:.2%}".format).rename(
            "% of LQ cells of cells with >250k reads"),
        place_holder.rename("line_break_3"),

        hq_metrics["hq_mean_quality"].apply(lambda x: round(x, 2)).rename("Mean quality of HQ cells"),
        hq_metrics["hq_median_quality"].apply(lambda x: round(x, 2)).rename("Median quality of HQ cells"),
        hq_metrics["hq_median_reads"].apply(lambda x: "{}k".format(int(x/1000))).rename("Median reads of HQ cells"),
        hq_metrics["hq_median_coverage_depth"].rename("Median coverge depth of HQ cells"),
        ]


    metrics = pd.concat(metrics_lst, axis=1)
    metrics = metrics.fillna(0)

    if "species" in df:
        fastqscreen_metrics_lst = [
        place_holder.rename("HQ cells (with quality >0.75)"),
        hq_metrics["hq_human_count"].rename("HQ human cells"),
        hq_metrics["hq_mouse_count"].rename("HQ mouse cells"),
        hq_metrics["hq_salmon_count"].rename("HQ salmon cells"),
        
        place_holder.rename("line_break_0"),
        place_holder.rename("Cells that have reads >250k"),
        non_dropout_metrics["human_count"].rename("human cells"),
        non_dropout_metrics["mouse_count"].rename("mouse cells"),
        non_dropout_metrics["salmon_count"].rename("salmon cells"),
        non_dropout_metrics["non_dropout_cells"].rename("Total cells that didn't drop out"),
        
        place_holder.rename("line_break_1"), 
        place_holder.rename("% of cells that have reads >250k"),
        non_dropout_metrics["human_ratio"].astype(float).map("{:.2%}".format).rename("% human cells"),
        non_dropout_metrics["mouse_ratio"].astype(float).map("{:.2%}".format).rename("% mouse cells"),
        non_dropout_metrics["salmon_ratio"].astype(float).map("{:.2%}".format).rename("% salmon cells"),
        place_holder.copy().apply(lambda x: "100%").rename("Percentage of non-dropout cells"),

        place_holder.rename("line_break_4"),

        hq_metrics["median_unmapped_ratio"].astype(float).map("{:.2%}".format).rename(
            "Median % of unmapped reads in HQ cells"),
        hq_metrics["median_unknown_reads_ratio"].astype(float).map("{:.2%}".format).rename(
            "Median % of unknown reads in HQ cells"),

        ]


        fastqscreen_metrics = pd.concat(fastqscreen_metrics_lst, axis=1)
        fastqscreen_metrics = fastqscreen_metrics.fillna(0)
    
        return metrics.T, fastqscreen_metrics.T
    else:
        return metrics.T, pd.DataFrame()

def generate_library_metrics(df, gc_data, reference_gc):
    cells_pass_df = df[df["quality"] >= 0.75]
    # Calculate total cells pass
    total_cells_pass = cells_pass_df.count()["cell_id"]

    if reference_gc is not None:
        gc_curve = get_gc_curve(gc_data)
        error_curve = reference_gc.subtract(gc_curve).abs()
        sum_of_diffs = float(error_curve.sum())
    else:
        sum_of_diffs = float('nan')

    df = pd.DataFrame.from_dict(
        {
            'Cell Yield': int(total_cells_pass),
            'sum of GC bias differences': sum_of_diffs,
        },
        orient='index'
    )
    df = df.fillna(0)

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
    # Compare with reference gc curve
    gc_curve.index = list(gc_curve.index.astype(int))

    if reference_gc is not None:
        reference_gc.index = list(reference_gc.index)
        plt.plot(reference_gc, color='r', label='reference curve')

    plt.grid(b=True, which='both', axis='both')
    plt.plot(gc_curve, 'o', mfc='none', label='Normalized Coverage')
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
        dict(selector="th", props=[('font-size', '12pt'), ("font-family", "Helvetica"),
                                   ("color", 'black'),
                                   ('background-color', 'rgb(232, 232, 232)'),
                                   ("text-align", "right"),
                                   ]),

        dict(selector="td", props=[("font-family", "Helvetica"),
                                   ("text-align", "right"),
                                   ('font-size', '11pt')
                                  ]),

        dict(selector="tr", props=[("line-height", "14px")]),
        dict(selector="caption", props=[("caption-side", "bottom")]),
    ]


def generate_html(dataframes, pngs, html_file):
    html_elements = []

    for header, df in dataframes:
        html_elements.append("<h3>{}</h3>\n".format(header))
        html_string = df.style.set_table_styles(pretty()).render()
        for place_holder in ["place_holder", "line_break_0", "line_break_1", "line_break_2", "line_break_3", "line_break_4"]:
            html_string = html_string.replace(place_holder, "<br>")
        html_elements.append(html_string)
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

    if reference_gc:
        reference_gc = load_reference(reference_gc)

    data = load_data(metrics)
    gc_data = load_data(gc_metrics, gc=True)

    qc_df, fastqscreen_df = generate_qc_table(data)
    plot_gc_curve(gc_data, reference_gc, gc_plot)
    plot_heatmap(data, heatmap)

    if not fastqscreen_df.empty:
        generate_html(
        [
            ('Metrics', qc_df),
            ('', fastqscreen_df),
        ],
        [
            ('Heatmap', heatmap),
            ('GC plot', gc_plot)
        ],
        html
    )
    else:
        generate_html(
        [
            ('Metrics', qc_df),
        ],
        [
            ('Heatmap', heatmap),
            ('GC plot', gc_plot)
        ],
        html
    )
