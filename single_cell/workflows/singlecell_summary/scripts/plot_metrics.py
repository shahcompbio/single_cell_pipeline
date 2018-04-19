'''
Plot sequencing metrics based on metric table and SampleSheet files.
'''
from __future__ import division

import argparse
import matplotlib
matplotlib.use('Agg')  # required for running on the cluster
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import warnings
import re

from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
matplotlib.rcParams['pdf.fonttype'] = 42


def parse_args():
    #=========================================================================
    # Read Command Line Input
    #=========================================================================
    parser = argparse.ArgumentParser()
    
    parser.add_argument('metric_table',
                        help='''Path to metric table for the run.''')
    
    parser.add_argument('out_file',
                        help='Path to output file where .pdf'
                        ' plots will be written.')
    
    parser.add_argument('--plot_title',
                        help='plot title to differentiate'
                        ' QC runs from full runs.')
    
    args = parser.parse_args()

    return args
#=========================================================================
# Functions
#=========================================================================

class PlotMetrics(object):
    def __init__(self, metrics, output, plot_title):
        self.metrics = metrics
        self.output = output
        self.plot_title = plot_title

    def add_barplot_labels(self, ax, labels, text_spacing, font_size):
        patch_width = list(set([x.get_width() for x in ax.patches]))[0]
    
        for plot_patch, label in zip(ax.patches, labels):
            patch_height = plot_patch.get_height()
    
            # if the bottom of rectangle is below 0
            if plot_patch.get_y() < 0:
                patch_height = -patch_height
    
            if np.isnan(patch_height):
                patch_height = 0
    
            if patch_height < 0:
                patch_height -= text_spacing
            elif patch_height > 0:
                patch_height += text_spacing
    
            ax.text(plot_patch.get_x() + patch_width / 2,
                    patch_height,
                    label,
                    ha='center',
                    va='bottom',
                    size=font_size
                    )
    
        return ax
    
    def plot_metric(self, df, metric, ylab, text_spacing, pdf, plot_title):
        sns.set(context='talk',
                style='ticks',
                font='Helvetica',
                rc={'axes.titlesize': 12,
                    'axes.labelsize': 12,
                    'xtick.labelsize': 12,
                    'ytick.labelsize': 12,
                    'legend.fontsize': 12})
    
        if text_spacing > max(df[metric]):
            warnings.warn('default text spacing is very high, overriding')
            text_spacing = 0.2 * max(df[metric])
    
        fig = plt.figure(figsize=(len(df['cell_id']) / 4, 5))
    
        ax = fig.gca()
    
        col = '#595959'
    
        sns.barplot(df['cell_id'], df[metric], color=col, ax=ax)
    
        column_labels = [str(x) for x in df['cell_call']]
    
        ax = self.add_barplot_labels(ax, column_labels, text_spacing, 12)
    
        ax.set_xlabel('Sample')
        ax.set_ylabel(ylab)
        sns.despine(offset=10, trim=True)
    
        sample_condition = [' (' + str(x) + ')' for x in
                            df['experimental_condition']]
        sample_labels = [x + y for x, y in zip(df['cell_id'], sample_condition)]
    
        ax.set_xticklabels(sample_labels)
    
        if plot_title:
            ax.set_title(plot_title, y=1.08, fontsize=10)
        plt.xticks(rotation=90)
    
        pdf.savefig(bbox_inches='tight', pad_inches=0.4)
        plt.close()
    
    def plot_metric_heatmap(self, df, metric, title, pdf, plot_title, size=72):
    
        # don't plot if we don't have the chip plate info
        # will only happen if merge_pipeline and plate info is not provided
        if df['sample_plate'].unique()[0] == ['R1-C1'] and \
                len(df['sample_plate'].unique()) == 1:
            return

        # add row and column cols to df
        df[['row', 'col']] = df.sample_plate.str.extract('[R](\d*)_[C](\d*)').apply(pd.Series).astype(int)

        #set size based on the R-C in sample_plate
        size = max(size, max(df.row), max(df.col))

        matrix = np.empty((size, size,))
        matrix[:] = np.nan

        sns.set(context='talk',
                style='darkgrid',
                font='Helvetica',
                rc={'axes.titlesize': 9,
                    'axes.labelsize': 6,
                    'xtick.labelsize': 6,
                    'ytick.labelsize': 6,
                    'legend.fontsize': 6})
    
        _ = plt.figure(figsize=(15, 12))
    
        tick_labels = [x + 1 for x in range(size)]

        well_labels=False
        #if the cell call matches the C[1-9]+ format then add cell calls to the plot
        if all([re.match("C[0-9]+$|NTC", cc) for cc in df['cell_call']]):
            well_labels = np.empty((size, size,))
            well_labels[:] = np.nan
    
            for i in range(len(df)):
                row_idx = int(df.ix[i, 'sample_plate'].split('_')[0].replace('R', ''))
                col_idx = int(df.ix[i, 'sample_plate'].split('_')[1].replace('C', ''))

                label_value = df.ix[i]["cell_call"]
                label_value = 0 if label_value == 'NTC' else int(label_value.replace('C', ''))
                well_labels[row_idx - 1, col_idx - 1] = label_value
    
        for i in range(len(df)):
            row_idx = int(df.ix[i, 'sample_plate'].split('_')[0].replace('R', ''))
            col_idx = int(df.ix[i, 'sample_plate'].split('_')[1].replace('C', ''))
            matrix_value = float(df.ix[i, metric])
            matrix[row_idx - 1, col_idx - 1] = matrix_value

       

        #raise Exception(set([c for v in np.isnan(matrix) for c in v]))
        try:
            sns.heatmap(matrix,
                        xticklabels=tick_labels,
                        yticklabels=tick_labels,
                        linewidths=0.6,
                        square=True,
                        cbar=True,
                        annot=well_labels,
                        annot_kws={'size': 6})
        
            plt.title(title + '(' + plot_title + ')')
        
            pdf.savefig(bbox_inches='tight', pad_inches=0.4)
        except ValueError:
            warnings.warn("Couldn't generate plot")
    
        plt.close()
    
    def plot_metric_factorplot(self, df, metric, ylab, pdf, plot_title):
        df_melt = pd.melt(df, id_vars=['cell_id', 'experimental_condition',
                                       'cell_call'], value_vars=[metric])
    
        expt_conditions_ordered = sorted(df['experimental_condition'].unique())
        cell_calls_ordered = sorted(df['cell_call'].unique())
    
        num_cell_calls = len(cell_calls_ordered)
        tableau_10_medium = ['#729ece', '#ff9e4a', '#67bf5c', '#ed665d', '#ad8bc9',
                             '#a8786e', '#ed97ca', '#a2a2a2', '#cdcc5d', '#6dccda']
        cols = tableau_10_medium[0:num_cell_calls]
    
        num_libs = []
    
        for condition in expt_conditions_ordered:
            num_call = []
    
            for call in cell_calls_ordered:
                num_call.append(str(
                    len(df[(df['experimental_condition'] == condition) & (df['cell_call'] == call)])))
    
            num_call = '\n'.join([','.join(num_call[x:x + 4])
                                  for x in range(0, len(num_call), 4)])
            num_libs.append(num_call)
    
        condition_labels = [x + '\n(n=' + y + ')'
                            for x, y in zip(expt_conditions_ordered, num_libs)]
    
    
        #need to scale fontsize if num_ec is low, or else tight_layout fails
        fontsize = min(12, 6 * len(expt_conditions_ordered))
        sns.set(context='talk',
                style='ticks',
                font='Helvetica',
                rc={'axes.titlesize': fontsize,
                    'axes.labelsize': fontsize,
                    'xtick.labelsize': fontsize,
                    'ytick.labelsize': fontsize,
                    'legend.fontsize': fontsize})
    
        fig_height = 6
        fig_width = len(expt_conditions_ordered) * 1.5
    
        aspect = max(0.5,  fig_width / fig_height)
    
        fig = plt.figure(figsize=(fig_width, fig_height))
        _ = fig.gca()
    
        g = sns.factorplot('experimental_condition',
                           'value',
                           'cell_call',
                           df_melt,
                           kind='box',
                           order=expt_conditions_ordered,
                           hue_order=cell_calls_ordered,
                           palette=cols,
                           legend=False,
                           size=fig_height,
                           aspect=aspect)
    
        plt.legend()
    
        g.set_axis_labels('Experimental condition', ylab)
        g.set_xticklabels(condition_labels, rotation=60)
    
        if plot_title:
            g.fig.suptitle(plot_title, fontsize=10)
    
        sns.despine(trim=True)
    
        pdf.savefig(bbox_inches='tight', pad_inches=0.4)
        plt.close()
    
    def sort_samples(self, df):
    
        # sort the df, by row and then by col
    
    
        df['row'] = df['cell_id'].str.extract('.*-R([0-9]*)-C[0-9]*').astype(int)
        df['col'] = df['cell_id'].str.extract('.*-R[0-9]*-C([0-9]*)').astype(int)
        df = df.sort_index(by=['row', 'col'], ascending=[True, True])
    
        return df
    
    #=========================================================================
    # Run script
    #=========================================================================
    
    def main(self, ):
        df = pd.read_csv(self.metrics)
    
        df = self.sort_samples(df)
    
        with PdfPages(self.output) as pdf:
            self.plot_metric(df, 'log_likelihood', 'Log Likelihood', 500,
                        pdf, self.plot_title,)
            self.plot_metric(df, 'mad_neutral_state', 'Mad Neutral State', 0.01,
                        pdf, self.plot_title,)
            self.plot_metric(df, 'MSRSI_non_integerness', 'MSRSI Non Integerness', 0.01,
                        pdf, self.plot_title,)
            self.plot_metric(df, 'MBRSI_dispersion_non_integerness', 'MBRSI Dispersion Non Integerness', 0.01,
                        pdf, self.plot_title,)
            self.plot_metric(df, 'MBRSM_dispersion', 'MBRSM Dispersion', 0.01,
                        pdf, self.plot_title,)

            self.plot_metric_heatmap(df, 'log_likelihood', 'Log Likelihood',
                                pdf, self.plot_title)
            self.plot_metric_heatmap(df, 'mad_neutral_state', 'Mad Neutral State',
                                pdf, self.plot_title)
            self.plot_metric_heatmap(df, 'MSRSI_non_integerness', 'MSRSI Non Integerness',
                                pdf, self.plot_title)
            self.plot_metric_heatmap(df, 'MBRSI_dispersion_non_integerness', 'MBRSI Dispersion Non Integerness',
                                pdf, self.plot_title)
            self.plot_metric_heatmap(df, 'MBRSM_dispersion', 'MBRSM Dispersion',
                                pdf, self.plot_title)

            self.plot_metric_factorplot(df, 'log_likelihood', 'Log Likelihood',
                                   pdf, self.plot_title)
            self.plot_metric_factorplot(df, 'mad_neutral_state', 'Mad Neutral State',
                                   pdf, self.plot_title)
            self.plot_metric_factorplot(df, 'MSRSI_non_integerness', 'MSRSI Non Integerness',
                                   pdf, self.plot_title)
            self.plot_metric_factorplot(df, 'MBRSI_dispersion_non_integerness', 'MBRSI Dispersion Non Integerness',
                                   pdf, self.plot_title)
            self.plot_metric_factorplot(df, 'MBRSM_dispersion', 'MBRSM Dispersion',
                                   pdf, self.plot_title)

if __name__ == '__main__':
    args = parse_args()
    
    plotter = PlotMetrics(args.metric_table, args.out_file, args.plot_title)
    plotter.main()
