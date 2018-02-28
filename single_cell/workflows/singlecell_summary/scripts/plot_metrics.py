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
    
    parser.add_argument('--gcbias_matrix',
                        help='gcbias matrix file')
    
    parser.add_argument('--gc_content_data',
                        help='gcbias matrix file')
    
    
    args = parser.parse_args()

    return args
#=========================================================================
# Functions
#=========================================================================

class PlotMetrics(object):
    def __init__(self, metrics, output, plot_title, gcbias_matrix, gc_content):
        self.metrics = metrics
        self.output = output
        self.plot_title = plot_title
        self.gcbias_matrix= gcbias_matrix
        self.gc_content = gc_content

    def add_legend(self, ax, labels, colours, num_columns, typ='rectangle',
                   location='upper center'):
        object_list = []
        for col in colours:
            if typ == 'rectangle':
                object_list.append(Rectangle((0, 0), 1, 1, facecolor=col,
                                             edgecolor='none'))
    
            elif typ == 'circle':
                object_list.append(Line2D(range(1), range(1), color=col,
                                          marker='o', markersize=15,
                                          linewidth=0)
                                   )
    
            else:
                warnings.warn('Legend type must be one of: rectangle, circle.')
    
        return ax.legend(tuple(object_list), tuple(labels), loc=location,
                         ncol=num_columns)
    
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
    
    def plot_metric_fraction_total(self, df, metric, metric_label, pdf,
                                   plot_title, from_top=False):
        sns.set(context='talk',
                style='ticks',
                font='Helvetica',
                rc={'axes.titlesize': 12,
                    'axes.labelsize': 12,
                    'xtick.labelsize': 12,
                    'ytick.labelsize': 12,
                    'legend.fontsize': 12})
    
        fig = plt.figure(figsize=(len(df['cell_id']) / 4, 5))
    
        ax = fig.gca()
    
        total = df['total_reads'] / 1000000
        fraction = df[metric] / 1000000
    
        col_total = '#cfcfcf'
        col_fraction = '#595959'
    
        if from_top:
            fraction = total - fraction
            sns.barplot(df['cell_id'], total, color=col_fraction, ax=ax)
            sns.barplot(df['cell_id'], fraction, color=col_total, ax=ax)
        else:
            sns.barplot(df['cell_id'], total, color=col_total, ax=ax)
            sns.barplot(df['cell_id'], fraction, color=col_fraction, ax=ax)
    
        column_labels = [str(x) for x in df['cell_call']]
    
        ax = self.add_barplot_labels(ax, column_labels, 0.05, 12)
    
        ax.set_xlabel('Sample')
        ax.set_ylabel('Number of reads (millions)')
        sns.despine(offset=10, trim=True)
    
        sample_condition = [' (' + str(x) + ')' for x in
                            df['experimental_condition']]
        sample_labels = [x + y for x, y in zip(df['cell_id'], sample_condition)]
    
        ax.set_xticklabels(sample_labels)
    
        plt.xticks(rotation=90)
    
        if plot_title:
            ax.set_title(plot_title, y=1.08, fontsize=10)
    
        self.add_legend(ax, ['Total', metric_label], [col_total, col_fraction], 1,
                   location='upper right')
    
        pdf.savefig(bbox_inches='tight', pad_inches=0.4)
        plt.close()
    
    def plot_metric_fraction(self, df, numerator_metric, denominator_metric,
                             ylab, pdf, plot_title):
        sns.set(context='talk',
                style='ticks',
                font='Helvetica',
                rc={'axes.titlesize': 12,
                    'axes.labelsize': 12,
                    'xtick.labelsize': 12,
                    'ytick.labelsize': 12,
                    'legend.fontsize': 12})
    
        fig = plt.figure(figsize=(len(df['cell_id']) / 4, 5))
    
        ax = fig.gca()
    
        fraction = df[numerator_metric] / df[denominator_metric]
    
        col_fraction = '#595959'
    
        sns.barplot(df['cell_id'], fraction, color=col_fraction, ax=ax)
    
        column_labels = [str(x) for x in df['cell_call']]
    
        ax = self.add_barplot_labels(ax, column_labels, 0.015, 12)
    
        ax.set_xlabel('Sample')
        ax.set_ylabel(ylab)
        plt.ylim(0, 1)
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
    
    def get_cmap(self, metrics):
        """
        generate dict with all cell calls and a randomly assigned
        color for each
        """
        ccs = set(metrics["cell_call"])
        cmap = sns.color_palette("deep", len(ccs))
    
        cmap = {cc: cm for cc, cm in zip(ccs, cmap)}
    
        return cmap
    
    
    def get_alpha(self, num_samples):
    
        alpha = 15.0/num_samples
    
        if alpha < 0.00001:
            return 0.0005
        elif alpha > 0.1:
            return 0.05
    
        return alpha
    
    
    
    def plot_gcbias_all(self, infile, pdf, plot_title, metrics, gcdata):
        df = pd.read_csv(infile)
        df = df.set_index("cell_id")
    
        cmap = self.get_cmap(metrics)

        samples = df.index
    
        plt.figure(figsize=(12,12))
        sns.set(context='talk',
                style='whitegrid',
                font='Helvetica',
                rc={'axes.titlesize': 9,
                    'axes.labelsize': 6,
                    'xtick.labelsize': 6,
                    'ytick.labelsize': 6,
                    'legend.fontsize': 6})
    
    
    
        alpha = self.get_alpha(len(samples))
    
        for samp in samples:
            cc = metrics[metrics['cell_id'] == samp]["cell_call"].iloc[0]
            plt.plot(range(0,101), df.loc[samp][map(str, range(0,101))], color=cmap[cc], alpha=alpha)
    
        if self.gc_content:
            ax = sns.barplot(x='gc', y='windows', data=gcdata,
                             color='#E7B591', ci=None)
            plt.setp(ax.patches, linewidth=0)
    
        ticks = np.arange(0, 100, 10)
        plt.xticks(ticks, map(str, ticks))
    
        plt.xlabel('GC% of 100 base windows')
        plt.ylabel('Normalized Coverage')
    
        plt.ylim((0, 2))
    
        # Plot the legend
        patches = [matplotlib.patches.Patch(color='#E7B591',
                                            label="Windows at GC%")]
        legend1 = plt.legend(handles=patches, bbox_to_anchor=(0, 0, 0.15, -0.15))
        plt.gca().add_artist(legend1)
    
        patches = [matplotlib.patches.Patch(color=v, label=k)
                   for k, v in cmap.iteritems()]
        plt.legend(handles=patches, bbox_to_anchor=(0, 0, 0.6, -0.15), ncol=6)
    
        pdf.savefig(bbox_inches='tight', pad_inches=0.4)
        plt.close()
    
    def plot_gcbias_by_ec(self, infile, pdf, plot_title, metrics, gcdata):
        """
        generate gcbias curves for all samples by ec and cc in legend
        """
        def get_samples_by_ec(metrics):
            """
            returns samples that belong to each ec
            """
            outdata = {}
            metrics = metrics.groupby("experimental_condition")
            for gc in metrics.groups.keys():
                vals = metrics.get_group(gc)["cell_id"]
    
                outdata[gc] = vals
            return outdata
    
        plt.figure(figsize=(12,12))
        sns.set(context='talk',
                style='whitegrid',
                font='Helvetica',
                rc={'axes.titlesize': 9,
                    'axes.labelsize': 6,
                    'xtick.labelsize': 6,
                    'ytick.labelsize': 6,
                    'legend.fontsize': 6})
    
        df = pd.read_csv(infile)
        df = df.set_index("cell_id")
    
        samps = get_samples_by_ec(metrics)
        cmap = self.get_cmap(metrics)
    
        #we dont want different alpha on different pages, so calculate it using max
        alpha = self.get_alpha(max([len(v) for v in samps.itervalues()]))
        for ec, samps in samps.iteritems():
            for samp in samps:
                cc = metrics[metrics['cell_id'] == samp]["cell_call"].iloc[0]
                plt.plot(range(0,101), df.loc[samp][map(str, range(0,101))], color=cmap[cc], alpha=alpha)
    
            if self.gc_content:
                ax = sns.barplot(x='gc', y='windows', data=gcdata,
                                 color='#E7B591', ci=None)
                plt.setp(ax.patches, linewidth=0)
    
            # Plot the legend
            patches = [matplotlib.patches.Patch(color='#E7B591',
                                                label="Windows at GC%")]
            legend1 = plt.legend(handles=patches,
                                 bbox_to_anchor=(0, 0, 0.15, -0.15))
            plt.gca().add_artist(legend1)
    
            patches = [matplotlib.patches.Patch(color=v, label=k)
                       for k, v in cmap.iteritems()]
            plt.legend(handles=patches, bbox_to_anchor=(0, 0, 0.6, -0.15), ncol=6)
    
            plt.xlabel('GC% of 100 base windows')
            plt.ylabel('Normalized Coverage')
            plt.title('experimental_condition:' + ec)
            ticks = np.arange(0, 100, 10)
            plt.xticks(ticks, map(str, ticks))
    
            plt.ylim((0, 2))
    
            pdf.savefig(bbox_inches='tight', pad_inches=0.4)
            plt.close()
    
    def plot_gcbias_by_ec_cc(self, infile, pdf, plot_title, metrics, gcdata):
        """
        generate gcbias curves for all samples by cc and ec
        """
        def get_samples_by_ec_cc(metrics):
            """
            returns samples that belong to each ec and cc
            """
            outdata = {}
            metrics = metrics.groupby(["experimental_condition", "cell_call"])
            for gc in metrics.groups.keys():
                vals = metrics.get_group(gc)["cell_id"]
    
                outdata[gc] = vals
            return outdata
    
        sns.set(context='talk',
                style='whitegrid',
                font='Helvetica',
                rc={'axes.titlesize': 9,
                    'axes.labelsize': 6,
                    'xtick.labelsize': 6,
                    'ytick.labelsize': 6,
                    'legend.fontsize': 6})
    
        df = pd.read_csv(infile)
        df = df.set_index("cell_id")
    
        samps = get_samples_by_ec_cc(metrics)
    
        #we dont want different alpha on different pages, so calculate it using max
        alpha = self.get_alpha(max([len(v) for v in samps.itervalues()]))
        for ec, samps in samps.iteritems():
            plt.figure(figsize=(12,12))
            for samp in samps:
                plt.plot(range(0,101), df.loc[samp][map(str, range(0,101))], color='#2098AE', alpha=alpha)
    
            if self.gc_content:
                ax = sns.barplot(x='gc', y='windows', data=gcdata,
                                 color='#E7B591', ci=None)
                plt.setp(ax.patches, linewidth=0)
    
            patches = [matplotlib.patches.Patch(color='#E7B591',
                                                label="Windows at GC%")]
            plt.legend(handles=patches, bbox_to_anchor=(0, 0, 0.5, -0.15))
    
            plt.ylim((0, 2))
            ticks = np.arange(0, 100, 10)
            plt.xticks(ticks, map(str, ticks))
            plt.xlabel('GC% of 100 base windows')
            plt.ylabel('Normalized Coverage')
            plt.title('experimental_condition: %s Cell Call %s' % (ec[0], ec[1]))
            pdf.savefig(bbox_inches='tight', pad_inches=0.4)
            plt.close()
    
    def read_gc_content(self, ):
        """
        read the file with the windows at %GC content information
        scaling is based on the Picard tools scaling method
        see https://github.com/broadinstitute/picard/blob/
        3dcadca7bf38a8cc9f6922b2334d082875899766/src/main/resources/picard/analysis/gcBias.R
        """
        if not self.gc_content:
            return None
    
        data = pd.read_csv(self.gc_content)
    
        win_ratio = 0.5 / max(data['windows'])
    
        data['windows'] = data['windows'] * win_ratio
    
        return data
    
    def plot_by_barcodes(self, df, metric, ylab, xlab, pdf, plot_title):
    
        sns.set(context='talk',
                style='ticks',
                font='Helvetica',
                rc={'axes.titlesize': 12,
                    'axes.labelsize': 12,
                    'xtick.labelsize': 12,
                    'ytick.labelsize': 12,
                    'legend.fontsize': 12})
    
        df_melt = pd.melt(df, id_vars=['cell_id', 'i5_barcode', 'i7_barcode'],
                          value_vars=[metric],)
    
        fig_height = 6
        fig_width = len(set(df_melt[xlab]))
    
        fig = plt.figure(figsize=(fig_width, fig_height))
        _ = fig.gca()
    
        g = sns.boxplot(x=df_melt[xlab], y=df_melt['value'],
                        orient='v', color='gray')
    
        g.set_ylabel(ylab)
        g.set_xlabel(xlab)
    
        if plot_title:
            g.set_title(plot_title, fontsize=10)
    
        plt.legend()
    
        sns.despine(trim=True)
    
        pdf.savefig(bbox_inches='tight', pad_inches=0.4)
        plt.close()
    
    #=========================================================================
    # Run script
    #=========================================================================
    
    def main(self, ):
        df = pd.read_csv(self.metrics)
    
        df = self.sort_samples(df)
    
        with PdfPages(self.output) as pdf:
            self.plot_metric_fraction_total(df, 'total_mapped_reads', 'Mapped', pdf,
                                       self.plot_title, from_top=False)
            self.plot_metric_fraction_total(df, 'total_duplicate_reads', 'Duplicates', pdf,
                                       self.plot_title, from_top=True)
            self.plot_metric_fraction_total(df, 'total_properly_paired', 'Properly paired', pdf,
                                       self.plot_title, from_top=False)

            self.plot_metric_fraction(df, 'total_mapped_reads', 'total_reads', 'Fraction mapped of total',
                                 pdf, self.plot_title,)
            self.plot_metric_fraction(df, 'total_duplicate_reads', 'total_mapped_reads', 'Fraction duplicates of mapped',
                                 pdf, self.plot_title,)
            self.plot_metric_fraction(df, 'total_properly_paired', 'total_mapped_reads',
                                 'Fraction properly paired of mapped', pdf, self.plot_title,)

            self.plot_metric(df, 'coverage_depth', 'Coverage depth', 0.0015,
                        pdf, self.plot_title,)
            self.plot_metric(df, 'coverage_breadth', 'Coverage breadth', 0.0015,
                        pdf, self.plot_title,)
            self.plot_metric(df, 'mean_insert_size', 'Mean insert size', 0.05,
                        pdf, self.plot_title,)
            self.plot_metric(df, 'median_insert_size', 'Median insert size', 0.05,
                        pdf, self.plot_title,)
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

            self.plot_metric_heatmap(df, 'total_reads', 'Total reads',
                                pdf, self.plot_title)
            self.plot_metric_heatmap(df, 'percent_duplicate_reads', 'Percent duplicate reads',
                                pdf, self.plot_title)
            self.plot_metric_heatmap(df, 'coverage_depth', 'Coverage depth',
                                pdf, self.plot_title)
            self.plot_metric_heatmap(df, 'coverage_breadth', 'Coverage breadth',
                                pdf, self.plot_title)
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

            self.plot_metric_factorplot(df, 'total_reads', 'Total reads',
                                   pdf, self.plot_title)
            self.plot_metric_factorplot(df, 'total_mapped_reads', 'Total mapped reads',
                                   pdf, self.plot_title)
            self.plot_metric_factorplot(df, 'percent_duplicate_reads',
                                   'Percent duplicate reads', pdf, self.plot_title)
            self.plot_metric_factorplot(df, 'coverage_depth', 'Coverage depth',
                                   pdf, self.plot_title)
            self.plot_metric_factorplot(df, 'coverage_breadth', 'Coverage breadth',
                                   pdf, self.plot_title)
            self.plot_metric_factorplot(df, 'mean_insert_size', 'Mean insert size',
                                   pdf, self.plot_title)
            self.plot_metric_factorplot(df, 'median_insert_size', 'Median insert size',
                                   pdf, self.plot_title)
            self.plot_metric_factorplot(df, 'estimated_library_size',
                                   'Picard estimated library size', pdf, self.plot_title)
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
            self.plot_by_barcodes(df, 'total_reads', 'Total Reads', 'i5_barcode',
                             pdf, self.plot_title)
            self.plot_by_barcodes(df, 'total_reads', 'Total Reads', 'i7_barcode',
                             pdf, self.plot_title)

            if self.gcbias_matrix:
                gcdata = self.read_gc_content()
                self.plot_gcbias_all(self.gcbias_matrix, pdf, self.plot_title,
                                df, gcdata)
                self.plot_gcbias_by_ec(self.gcbias_matrix, pdf, self.plot_title,
                                  df, gcdata)
                self.plot_gcbias_by_ec_cc(self.gcbias_matrix, pdf, self.plot_title,
                                     df, gcdata)

if __name__ == '__main__':
    args = parse_args()
    
    plotter = PlotMetrics(args.metric_table, args.out_file, args.plot_title, args.gcbias_matrix, args.gc_content_data)
    plotter.main()
