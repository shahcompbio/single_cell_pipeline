'''
Created on Jul 24, 2017

@author: dgrewal
'''
import pandas as pd
from scripts import SummaryMetrics
from scripts import PlotKernelDensity
from scripts import MergeFiles
from scripts import PlotHeatmap
from scripts import PlotMetrics
from scripts import PlotPcolor


def concatenate_csv(in_filenames, out_filename, nan_val='NA'):
    data = []
    for _, in_filename in in_filenames.iteritems():
        with open(in_filename) as f:
            first_line = f.readline()
            if len(first_line) == 0:
                continue
        data.append(pd.read_csv(in_filename, dtype=str))
    data = pd.concat(data, ignore_index=True)
    data = data.fillna(nan_val)
    data.to_csv(out_filename, index=False)


def merge_csv(in_filenames, out_filename, how, on, nan_val='NA'):
    data = []
    for _, in_filename in in_filenames.iteritems():
        with open(in_filename) as f:
            first_line = f.readline()
            if len(first_line) == 0:
                continue
        data.append(pd.read_csv(in_filename, dtype=str))

    data = merge_frames(data, how, on)
    data = data.fillna(nan_val)
    data.to_csv(out_filename, index=False)


def merge_frames(frames, how, on):
    '''
    annotates input_df using ref_df
    '''

    if ',' in on:
        on = on.split(',')

    if len(frames) == 1:
        return frames[0]
    else:
        left = frames[0]
        right = frames[1]
        merged_frame = pd.merge(left, right,
                                how=how,
                                on=on)
        for frame in frames[2:]:
            merged_frame = pd.merge(merged_frame, frame,
                                    how=how,
                                    on=on)
        return merged_frame


def get_summary_metrics(infile, output):
    summ = SummaryMetrics(infile, output)
    summ.main()


def plot_kernel_density(infile, output, sep, colname, plot_title):
    plot = PlotKernelDensity(infile, output, sep, colname, plot_title)
    plot.main()

def merge_tables(infile, output, typ, sep,
                 merge_type, key_cols, nan_val):

    m = MergeFiles(infile, output, typ, sep,
                   merge_type, key_cols, nan_val)
    m.main()


def plot_metrics(metrics, output, plot_title, gcbias_matrix, gc_content):
    plot = PlotMetrics(metrics, output, plot_title, gcbias_matrix, gc_content)
    plot.main()


def plot_heatmap(infile, metrics, order_data, output, plot_title=None,
                 colname=None, plot_by_col=None, numreads_threshold=None,
                 mad_threshold=None, chromosomes=None):

    plot = PlotHeatmap(infile, metrics, order_data, output, plot_title=plot_title,
                       colname=colname, plot_by_col=plot_by_col,
                       numreads_threshold=numreads_threshold,
                       mad_threshold=mad_threshold, chromosomes=chromosomes)
    plot.main()



def plot_pcolor(infile, metrics, order_data, output, plot_title=None,
                 colname=None, plot_by_col=None, numreads_threshold=None,
                 mad_threshold=None, chromosomes=None):

    plot = PlotPcolor(infile, metrics, order_data, output, plot_title=plot_title,
                       colname=colname, plot_by_col=plot_by_col,
                       numreads_threshold=numreads_threshold,
                       mad_threshold=mad_threshold, chromosomes=chromosomes)
    plot.main()
