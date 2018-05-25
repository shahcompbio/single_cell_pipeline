'''
Created on Nov 16, 2016

@author: dgrewal
'''
from __future__ import division

import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from matplotlib.backends.backend_pdf import PdfPages
import utils as utl
import matplotlib.gridspec as gridspec

from matplotlib.colors import rgb2hex


import numpy as np

lowess = sm.nonparametric.lowess

matplotlib.rcParams['pdf.fonttype'] = 42


def parse_args():

    ann_cols = ['pick_met', 'condition', 'sample_type',
                'median_hmmcopy_reads_per_bin', 'mad_neutral_state',
                'MSRSI_non_integerness']

    #=========================================================================
    # Read Command Line Input
    #=========================================================================
    parser = argparse.ArgumentParser()

    parser.add_argument('--corrected_reads',
                        required=True,
                        help='''Path to HMMcopy corrected reads output .csv file.''')

    parser.add_argument('--segments',
                        required=True,
                        help='''Path to HMMcopy segments output .csv file.''')

    parser.add_argument('--params',
                        required=True,
                        help='''Path to HMMcopy params output .csv file.''')

    parser.add_argument('--quality_metrics',
                        required=True,
                        help='''Optional quality metrics file for the run, with 'mad_neutral_state' column.''')

    parser.add_argument('--ref_genome',
                        required=True,
                        help='''Path to reference genome used for alignment.''')

    parser.add_argument('--num_states', type=int, default=7,
                        help='''Number of states used to run HMMcopy, default 7.''')

    parser.add_argument('--bias_output',
                        required=True,
                        help='''Path to HMMcopy bias reads output .pdf file.''')

    parser.add_argument('--segs_output',
                        required=True,
                        help='''Path to HMMcopy segs reads output .pdf file.''')

    parser.add_argument('--sample_id',
                        required=True,
                        help='''title of the plots''')

    parser.add_argument('--multipliers',
                        required=True,
                        nargs="*",
                        help='''scales to plot''')

    parser.add_argument('--annotation_cols',
                        default=ann_cols,
                        help='''title of the plots''')

    args = parser.parse_args()
    return args


class GenHmmPlots(object):
    """
    generate the reads, bias and segment plots
    """

    def __init__(self, reads, segments, params, metrics,
                 ref_genome, segs_out, bias_out,
                 sample_id, multipliers, **kwargs):

        self.multipliers = multipliers
        self.reads = reads
        self.segments = segments
        self.params = params
        self.metrics = metrics
        self.ref_genome = ref_genome
        self.sample_id = sample_id

        self.num_states = kwargs.get('num_states')
        if not self.num_states:
            self.num_states = 7

        self.annotation_cols = kwargs.get('annotation_cols')

        if not self.annotation_cols:
            self.annotation_cols = ['pick_met', 'condition',
                                    'sample_type', 'coverage_depth',
                                    'mad_neutral_state', 'MSRSI_non_integerness']

        self.segs_pdf = segs_out
        self.bias_pdf = bias_out

        self.sample_info = kwargs.get("sample_info")

    def __enter__(self):
        self.reads_store = pd.HDFStore(self.reads, 'r')
        self.segments_store = pd.HDFStore(self.segments, 'r')
        self.metrics_store = pd.HDFStore(self.metrics, 'r')
        self.params_store = pd.HDFStore(self.params, 'r')

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.reads_store.close()
        self.segments_store.close()
        self.params_store.close()
        self.metrics_store.close()

    def read_metrics(self, cell_id, multiplier):
        """

        """

        return self.metrics_store[
            "/hmmcopy/metrics/{}/{}".format(cell_id, multiplier)]

    def read_params(self, cell_id, multiplier):
        """

        """
        return self.params_store[
            "/hmmcopy/params/{}/{}".format(cell_id, multiplier)]

    def read_corrected_reads(self, cell_id, multiplier):
        """

        """

        df = self.reads_store[
            "/hmmcopy/reads/{}/{}".format(cell_id, multiplier)]

        df = utl.normalize_reads(df)
        df = utl.compute_chromosome_coordinates(df, self.ref_genome)

        return df

    def read_segments(self, cell_id, multiplier):
        """

        """

        df = self.segments_store[
            "/hmmcopy/segments/{}/{}".format(cell_id, multiplier)]

        chromosomes = map(str, range(1, 23)) + ['X', 'Y']
        df["chr"] = pd.Categorical(df["chr"], chromosomes)
        df = df.sort_values(['chr', 'start', 'end'])

        df = utl.compute_chromosome_coordinates(df, self.ref_genome)

        return df

    def get_annotations(self, metrics):
        annotations = []
        for colname in self.annotation_cols:

            if colname in metrics:
                val = metrics[colname].iloc[0]
                if not isinstance(val, str):
                    val = '%.3f' % val
            elif self.sample_info and colname in self.sample_info:
                val = self.sample_info[colname]
            else:
                val = 'NA'

            annstr = "{} : {}".format(colname, val)
            annotations.append(annstr)
        return annotations

    def add_annotations(
            self, fig, annotations, pos=(0.5, 0.02), fontsize=None):
        annotations = ', '.join(annotations)

        fig.text(pos[0], pos[1], annotations, fontsize=fontsize,
                 horizontalalignment='center',
                 verticalalignment='center',)

    def plot_bias(self, pdfout):
        """
        """

        sns.set(context='talk',
                style='ticks',
                font='Helvetica',
                rc={'axes.titlesize': 12,
                    'axes.labelsize': 12,
                    'xtick.labelsize': 12,
                    'ytick.labelsize': 12,
                    'legend.fontsize': 12,
                    'font.size': 12})

        df = self.read_corrected_reads(self.sample_id, 0)

        metrics = self.read_metrics(self.sample_id, 0)
        annotations = self.get_annotations(metrics)

        df_ideal = df[df['ideal'] == True]

        col = '#006ba4'

        fig = plt.figure(figsize=(12, 15))

        ax1 = plt.subplot(311)
        ax1.set_ylabel('Read count')
        ax1.set_xlabel('GC Content')
        if df is not None:
            ax1.scatter(
                df_ideal['gc'],
                df_ideal['reads'],
                edgecolors=col,
                facecolors='none',
                alpha=0.1,
                rasterized=True)

            mapp = {
                gc: pred for gc,
                pred in zip(
                    df_ideal["gc"],
                    df_ideal["modal_curve"])}
            x = sorted(df_ideal["gc"])
            y = [mapp[v] for v in x]
            plt.plot(x, y)

        ax2 = plt.subplot(312, sharex=ax1)
        ax2.set_xlabel('GC content')
        ax2.set_ylabel('Normalized read count')
        not_null = df_ideal['cor_gc'].notnull()
        if df is not None:
            ax2.scatter(
                df_ideal['gc'][not_null],
                df_ideal['cor_gc'][not_null],
                edgecolors=col,
                facecolors='none',
                alpha=0.1,
                rasterized=True)

        ax3 = plt.subplot(313)
        ax3.set_xlabel('Mappability')
        ax3.set_ylabel('Read Count')
        if df is not None:
            ax3.scatter(
                df_ideal['map'],
                df_ideal['reads'],
                edgecolors=col,
                facecolors='none',
                alpha=0.1,
                rasterized=True)

        self.add_annotations(fig, annotations, fontsize=8)
        fig.suptitle(self.sample_id, fontsize=12)
        sns.despine(offset=10, trim=True)

        plt.tight_layout(rect=(0, 0.05, 1, 0.95))
        pdfout.savefig(fig, pad_inches=0.2)
        plt.close()

    def get_colors(self, num_states):

        color_reference = {0: '#3498DB', 1: '#85C1E9', 2: '#808080'}


        low_max = 3 + int((num_states-3)/2) + 1
        hi_max = num_states + 1
        low_states = np.arange(3, low_max)
        hi_states = np.arange(low_max, hi_max)

        low_cmap = matplotlib.cm.get_cmap('OrRd', low_max+1)
        hi_cmap = matplotlib.cm.get_cmap('RdPu', hi_max+1)

        for cn_level in low_states:
            rgb = low_cmap(int(cn_level))[:3]
            color_reference[cn_level] = rgb2hex(rgb)

        for cn_level in hi_states:
            rgb = hi_cmap(int(cn_level))[:3]
            color_reference[cn_level] = rgb2hex(rgb)

        labels = [str(x) for x in range(num_states + 1)]
        labels[-1] = labels[-1] + ' or more'

        return color_reference


    def plot_segments(self, pdfout, remove_y=False):

        num_plots = len(self.multipliers)

        sns.set(context='talk',
                style='ticks',
                font='Helvetica',
                rc={'axes.titlesize': 6 * num_plots,
                    'axes.labelsize': 6 * num_plots,
                    'xtick.labelsize': 6 * num_plots,
                    'ytick.labelsize': 6 * num_plots,
                    'legend.fontsize': 6 * num_plots,
                    'font.size': 6 * num_plots})

        linewidth = 0.5 * num_plots
        scatter_size = 4 * num_plots

        height_plot = 20 * num_plots
        width_plot = 10 * num_plots

        # standard: 15,4
        # SA501X3F xenograft heatmap: 20.4, 4
        fig = plt.figure(figsize=(height_plot, width_plot))

        gs = gridspec.GridSpec(num_plots, 3, width_ratios=[20, 3, 1], wspace=0)

        # add annotations to plot
        metrics = self.read_metrics(self.sample_id, 0)
        annotations = self.get_annotations(metrics)
        self.add_annotations(fig, annotations)

        for i, multiplier in enumerate(self.multipliers):

            reads = self.read_corrected_reads(self.sample_id, multiplier)
            segments = self.read_segments(self.sample_id, multiplier)
            params = self.read_params(self.sample_id, multiplier)
            metrics = self.read_metrics(self.sample_id, multiplier)

            if reads is not None and remove_y:
                reads = reads[reads['chr'] != 'Y']

            ax = plt.subplot(gs[i, 0])
            ax = utl.create_chromosome_plot_axes(ax, self.ref_genome)

            ax.set_title("multiplier: {}".format(multiplier))

            num_states = np.nanmax(reads.state)
            cmap = self.get_colors(num_states)

            # we pass None if we don't have data
            if reads is not None and segments is not None:
                cols = reads["state"].replace(cmap)
                cols = cols[~reads['copy'].isnull()]

                df = reads[np.isfinite(reads['copy'])]

                plt.scatter(
                    df['plot_coord'],
                    df['copy'],
                    facecolors=cols,
                    edgecolors='none',
                    s=scatter_size,
                    rasterized=True)

                x, y = utl.get_segment_start_end(segments, remove_y)
                plt.plot(x, y, color='black', linewidth=linewidth)

            ylim = np.nanpercentile(reads['copy'], 99)
            if not np.isfinite(ylim):
                ylim = ax.get_ylim()[1]
            ylim = int(ylim) + 1
            ylim = max(3, ylim)

            if ylim <=10:
                tick_step = 1
            elif ylim > 10 and ylim <= 30:
                tick_step = 2
            else:
                tick_step = 5

            #make ylim a multiple of tick_step
            ylim = int((ylim+tick_step)/tick_step)*tick_step

            ax.set_ylim((0, ylim))

            ax.set_yticks(np.arange(0, ylim, tick_step))

            sns.despine(offset=10, trim=True)
            ax.tick_params(axis='x', which='minor', pad=9.1)

            ax1 = plt.subplot(gs[i, 1], sharey=ax)
            ax1.set_ylim((0, ylim))

            sns.despine(offset=10)

            ax1.yaxis.set_visible(False)
            ax1.spines['left'].set_visible(False)

            self.plot_dist(
                ax1,
                df,
                params,
                cmap,
                num_states,
                vertical=True)

            ax1.set_xlim((0, 2.5))


            ax.set_ylabel('Copy number')

            if i+1 == len(self.multipliers):
                ax.set_xlabel('Chromosome')
                ax1.set_xlabel("density")

            ax = utl.add_open_grid_lines(ax)

            ax2 = plt.subplot(gs[i, 2])
            ax2.axis('off')
            fig.legend = utl.add_legend(
                ax2,
                cmap,
                2,
                type='rectangle',
                location='upper center')


        fig.suptitle(self.sample_id, x=0.5, y=0.97)
        plt.tight_layout(rect=(0, 0.05, 1, 0.95))

        pdfout.savefig(fig)
        plt.close()

    def plot_dist(self, fig, reads, params,
                  cmap, num_states=7, vertical=False):

        # no data to plot
        if reads['copy'].isnull().all():
            return

        scale = (reads["copy"] / reads["copy"])
        scale = scale[~reads['copy'].isnull()].unique()
        # account for floating point errors
        scale = scale[np.isfinite(scale)]
        if not scale.any():
            return
        assert np.nanvar(scale) < 0.0001
        scale = scale[0]

        for state in range(0, num_states):
            color = cmap[state]

            data = reads[reads["state"] == state]["copy"]

            if np.isnan(data).all():
                continue

            if not data.empty:
                sns.kdeplot(data, bw="scott", kernel="epa",
                            shade=True, linewidth=0.5,
                            label=state, legend=False,
                            color=color, vertical=vertical)

            x = np.arange(0, np.nanmax(np.array(reads["copy"])), 0.01)
            mu = params[
                (params["parameter"] == "mus") & (
                    params["state"] == state)]["final"].iloc[0]
            lmbda = params[
                (params["parameter"] == "lambdas") & (
                    params["state"] == state)]["final"].iloc[0]
            nu = params[
                (params["parameter"] == "nus") & (
                    params["state"] == state)]["final"].iloc[0]

            y = utl.t_dist_pdf(x, mu, lmbda, nu)
            x = x * scale

            if vertical:
                plt.plot(y, x, color, linewidth=0.5, linestyle='--')
                fig.set_xlim(0, 4)
            else:
                plt.plot(x, y, color, linewidth=0.5, linestyle='--')
                fig.set_ylim(0, 4)


    def main(self):
        """
        main
        """

        with PdfPages(self.segs_pdf) as seg_pdf_out:
            self.plot_segments(seg_pdf_out)

        with PdfPages(self.bias_pdf) as bias_pdf_out:
            self.plot_bias(bias_pdf_out)


if __name__ == '__main__':
    args = parse_args()

    with GenHmmPlots(
        args.corrected_reads, args.segments, args.params, args.quality_metrics,
        args.ref_genome, args.segs_output, args.bias_output, args.sample_id,
        args.multipliers, annotation_cols=args.annotation_cols,
        num_states=args.num_states) as plot:
        plot.main()
