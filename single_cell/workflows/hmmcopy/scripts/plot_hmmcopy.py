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
import math
import warnings
from matplotlib.patches import Patch
import statsmodels.nonparametric.api as smnp

from matplotlib.colors import rgb2hex


import numpy as np

lowess = sm.nonparametric.lowess

matplotlib.rcParams['pdf.fonttype'] = 42

sns.set(context='talk',
        style='ticks',
        font='Helvetica',
        rc={'axes.titlesize': 12,
            'axes.labelsize': 12,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'legend.fontsize': 12,
            'font.size':12})


def parse_args():

    ann_cols = ['cell_call','experimental_condition', 'sample_type',
                'coverage_depth', 'mad_neutral_state', 'MSRSI_non_integerness']

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

    parser.add_argument('--sample_info',
                        help='''sample information''')

    parser.add_argument('--ref_genome',
                        required=True,
                        help='''Path to reference genome used for alignment.''')

    parser.add_argument('--num_states', type=int, default=7,
                        help='''Number of states used to run HMMcopy, default 7.''')

    parser.add_argument('--mad_threshold', type=float, default=0,
                        help='''all cells that have low MAD won't be plotted''')

    parser.add_argument('--reads_output',
                        required=True,
                        help='''Path to HMMcopy corrected reads output .pdf file.''')

    parser.add_argument('--bias_output',
                        required=True,
                        help='''Path to HMMcopy bias reads output .pdf file.''')

    parser.add_argument('--segs_output',
                        required=True,
                        help='''Path to HMMcopy segs reads output .pdf file.''')

    parser.add_argument('--params_output',
                        required=True,
                        help='''Path to HMMcopy segs reads output .pdf file.''')


    parser.add_argument('--plot_title',
                        help='''title of the plots''')

    parser.add_argument('--sample_id',
                        required=True,
                        help='''title of the plots''')

    parser.add_argument('--annotation_cols',
                        default=ann_cols,
                        help='''title of the plots''')

    args = parser.parse_args()
    return args


class GenHmmPlots(object):
    """
    generate the reads, bias and segment plots
    """

    def __init__(self, reads, segments, params, metrics, sample_info,
                 ref_genome, reads_out, segs_out, bias_out, params_out,
                 sample_id, **kwargs):

        self.reads = reads
        self.segments = segments
        self.params = params
        self.metrics = metrics
        self.ref_genome = ref_genome
        self.sample_id = sample_id

        self.sample_info = sample_info

        self.num_states = kwargs.get('num_states')
        if not self.num_states:
            self.num_states = 7

        self.mad_threshold = kwargs.get('mad_threshold')
        if not self.mad_threshold:
            self.mad_threshold = 0

        self.plot_title = kwargs.get('plot_title')
        if not self.plot_title:
            self.plot_title = ''

        self.annotation_cols = kwargs.get('annotation_cols')

        if not self.annotation_cols:
            self.annotation_cols = ['cell_call','experimental_condition',
                                    'sample_type', 'coverage_depth',
                                    'mad_neutral_state', 'MSRSI_non_integerness']

        self.reads_pdf, self.segs_pdf, self.bias_pdf, self.params_pdf = self.get_pdf_handles(
            reads_out, bias_out, segs_out, params_out)

    def load_data_pandas(self, infile):
        """

        """
        data = pd.read_csv(infile,
                           sep=',')

        return data

    def read_quality_metrics(self):
        """

        """

        df = self.load_data_pandas(self.metrics)

        df = df[df['cell_id'] == self.sample_id]

        return df


    def read_params(self):
        """

        """

        df = self.load_data_pandas(self.params)

        df = df[df['cell_id'] == self.sample_id]

        return df

    def read_corrected_reads(self):
        """

        """

        df = self.load_data_pandas(self.reads)

        df = utl.normalize_reads(df)
        df = utl.compute_chromosome_coordinates(df, self.ref_genome)

        return df

    def read_segments(self):
        """

        """

        df = self.load_data_pandas(self.segments)
        df = utl.compute_chromosome_coordinates(df, self.ref_genome)

        return df

    def get_pdf_handles(self, reads_out, bias_out, segs_out, params_out):
        """

        """

        reads_pdf = PdfPages(reads_out)
        bias_pdf = PdfPages(bias_out)
        segs_pdf = PdfPages(segs_out)
        params_pdf = PdfPages(params_out)

        return reads_pdf, segs_pdf, bias_pdf, params_pdf

    def get_annotations(self, metrics):
        annotations = []
        for colname in self.annotation_cols:

            if colname in metrics:
                val = metrics[colname].iloc[0]
                if not isinstance(val, str):
                    val = '%.3f' % val
            else:
                val = 'NA'

            annstr = "{} : {}".format(colname, val)
            annotations.append(annstr)
        return annotations

    def add_annotations(self, fig, annotations):
        annotations = ', '.join(annotations)

        fig.text(0.05, 0.02, annotations)

    def get_plot_title(self, sample_id, metrics):
        """

        """
        title_str = sample_id + self.plot_title
        return title_str

    def get_mad_score(self, sample_id, metrics):
        """
        """
        mad = metrics['mad_neutral_state'].iloc[0]
        return mad

    def gen_reads_plot(self, df, ax, typ='norm'):

        col = '#595959'

        ax = utl.create_chromosome_plot_axes(ax, self.ref_genome)

        # only attempt to plot if data is available
        if df is not None:
            plt.scatter(df['plot_coord'], df[typ], color=col, s=4)
            plt_lowess = lowess(
                df[typ],
                df['plot_coord'],
                frac=0.01,
                return_sorted=False)
            plt.plot(
                df['plot_coord'],
                plt_lowess,
                color='black',
                linewidth=1.2)

        if typ == 'norm':
            ax.set_ylabel('Normalized reads per bin')
        elif typ == 'cor_gc':
            ax.set_ylabel('GC corrected reads per bin')
        elif typ == 'cor_map':
            ax.set_ylabel('GC and mappability \n corrected reads per bin')
        ax.tick_params(axis='x', which='minor', pad=9.1)
        ax = utl.add_open_grid_lines(ax)

    def plot_corrected_reads(self, df, sample_id, title, annotations):
        """

        """
        fig = plt.figure(figsize=(15, 12))

        plt.subplot(3, 1, 1)
        ax = fig.gca()
        self.gen_reads_plot(df, ax, typ='norm')
        ax.set_title(title)
        ax.set_xlabel('')

        plt.subplot(3, 1, 2)
        ax = fig.gca()
        self.gen_reads_plot(df, ax, typ='cor_gc')
        ax.set_xlabel('')

        plt.subplot(3, 1, 3)
        ax = fig.gca()
        self.gen_reads_plot(df, ax, typ='cor_map')

        sns.despine(offset=10, trim=True)
        plt.tight_layout(rect=(0, 0.05, 1, 1))

        self.add_annotations(fig, annotations)
        self.reads_pdf.savefig(fig, pad_inches=0.2)
        plt.close()

    def plot_bias(self, df, sample_id, title, annotations):
        """
        """
        df_ideal = df[df['ideal'] == True]

        col = '#006ba4'

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
            2, 2, sharex='col', sharey='row', figsize=(15, 15))

        ax1.set_ylabel('Read count')
        ax1.set_title('Uncorrected')
        if df is not None:
            ax1.scatter(
                df_ideal['gc'],
                df_ideal['reads'],
                edgecolors=col,
                facecolors='none',
                alpha=0.1)

        ax2.set_title('Uncorrected')
        if df is not None:
            ax2.scatter(
                df_ideal['map'],
                df_ideal['reads'],
                edgecolors=col,
                facecolors='none',
                alpha=0.1)

        ax3.set_xlabel('GC content')
        ax3.set_ylabel('Normalized read count')
        ax3.set_title('GC corrected')
        not_null = df_ideal['cor_gc'].notnull()
        if df is not None:
            ax3.scatter(
                df_ideal['gc'][not_null],
                df_ideal['cor_gc'][not_null],
                edgecolors=col,
                facecolors='none',
                alpha=0.1)

        ax4.set_xlabel('Mappability')
        ax4.set_title('GC and mappability corrected')
        not_null = df_ideal['cor_map'].notnull()
        if df is not None:
            ax4.scatter(
                df_ideal['map'][not_null],
                df_ideal['cor_map'][not_null],
                edgecolors=col,
                facecolors='none',
                alpha=0.1)

        self.add_annotations(fig, annotations)
        fig.suptitle(title, fontsize=10)
        sns.despine(offset=10, trim=True)

        plt.tight_layout(rect=(0, 0.05, 1, 0.95))
        self.bias_pdf.savefig(fig, pad_inches=0.2)
        plt.close()

    def plot_segments(self, df, segments, plot_title,
                      annotations, num_states=7, remove_y=False):
        if df is not None and remove_y:
            df = df[df['chr'] != 'Y']

        # standard: 15,4
        # SA501X3F xenograft heatmap: 20.4, 4
        fig = plt.figure(figsize=(15, 4))
        ax = fig.gca()

        ax = utl.create_chromosome_plot_axes(ax, self.ref_genome)
        ax.set_title(plot_title)
        ax.set_ylabel('Copy number')

        segment_states = range(1, num_states + 1)
        segment_labels = [str(x) for x in range(num_states)]
        segment_labels[-1] = segment_labels[-1] + ' or more'
        segment_colours = [
            '#006ba4',
            '#5f9ed1',
            '#ababab',
            '#ffbc79',
            '#ff800e',
            '#c85200',
            '#8c3900']

        # we pass None if we don't have data
        if df is not None and segments is not None:
            cols = df['state']
            cols = cols.replace(segment_states, segment_colours)
            cols = cols[~df['copy'].isnull()]

            plt.scatter(
                df['plot_coord'],
                df['integer_copy_scale'],
                facecolors=cols,
                edgecolors='none',
                s=4)

            x, y = utl.get_segment_start_end(segments, remove_y)
            plt.plot(x, y, color='black', linewidth=1)

        #ax.set_ylim((0, plt.ylim()[1]))
        ax.set_ylim((0, 14))

        sns.despine(offset=10, trim=True)
        ax.tick_params(axis='x', which='minor', pad=9.1)

        ax.legend = utl.add_legend(
            ax,
            segment_labels,
            segment_colours,
            num_states,
            type='rectangle',
            location='upper center')
        ax = utl.add_open_grid_lines(ax)

        self.add_annotations(fig, annotations)

        plt.tight_layout(rect=(0, 0.05, 1, 1))
        self.segs_pdf.savefig(fig, pad_inches=0.2)
        plt.close()

    def get_sample_data(self, df, sample, norm=False):
        """
        """
        if sample not in df.groups:
            return None

        df = df.get_group(sample)
        if norm:
            df = utl.normalize_reads(df)
        df = utl.compute_chromosome_coordinates(df, self.ref_genome)

        return df

    def check_mad_score(self, sample, metrics):
        """

        """
        mad = self.get_mad_score(sample, metrics)

        # if mad_threshold is set to nonzero.
        # zero is defaults and means mad_threshold is not set. so no filtering
        if self.mad_threshold:
            if math.isnan(mad):
                return False

            if mad > self.mad_threshold:
                return False
        return True

    def check_info_columns(self, metrics):

        info_cols = ['cell_call', 'experimental_condition', 'sample_type']

        if not all(x in metrics.columns.values for x in info_cols):
            if not self.sample_info:
                warnings.warn("missing cell information in metrics, "
                              "plot titles might show NA. "
                              "please provide sample_info file")
            else:
                sample_info = self.load_data_pandas(self.sample_info)
                metrics = pd.merge(metrics, sample_info, on=['cell_id'])

        return metrics

    def get_cmap(self, num_states):
        """
        generating a custom heatmap 2:gray 0: blue 2+: reds
        """
        cmap = matplotlib.cm.get_cmap('coolwarm', num_states)
        colors = {}

        for i in range(cmap.N):
            # will return rgba, we take only first 3 so we get rgb
            rgb = cmap(i)[:3]
            colors[i+1] = rgb2hex(rgb)

        return colors

    def add_legend(self, fig, cmap):
        lgnd_patches = [Patch(color=c, label=k)
                        for k,c in cmap.iteritems()]
        plt.legend(handles=lgnd_patches,
                   bbox_to_anchor=(1, 1))

    def plot_params(self, reads, params, plot_title,
                      annotations, num_states=7):

        fig = plt.figure(figsize=(10, 10))
        cmap = self.get_cmap(num_states)

        if params.empty:
            return

        ylim_max = 0

        for state in range(1, num_states+1):
            color = cmap[state]
            data = reads[reads["state"] == state]["copy"]

            if not data.empty:
                sns.kdeplot(data, bw="scott",  kernel="epa",
                            shade=True, linewidth=0.5,
                            label=state, legend=False,
                            color = color)

            x = np.arange(0, np.nanmax(np.array(reads["copy"])), 0.01)
            mu = params[(params["parameter"]=="mus") & (params["state"] == state)]["final"].iloc[0]
            lmbda = params[(params["parameter"]=="lambdas") & (params["state"] == state)]["final"].iloc[0]
            nu = params[(params["parameter"]=="nus") & (params["state"] == state)]["final"].iloc[0]

            y = utl.t_dist_pdf(x, mu, lmbda, nu)
            plt.plot(x,y, color)
            ylim_max = max([ylim_max, max(y)])


        ylim_max = max(plt.ylim()[1], ylim_max)
        plt.ylim((0, ylim_max))
        plt.ylabel("density")
        plt.suptitle(plot_title)

        self.add_legend(fig, cmap)

        self.add_annotations(fig, annotations)
        plt.tight_layout(rect=(0, 0.05, 1, 0.95))
        self.params_pdf.savefig(fig)
        plt.close()

    def main(self):
        """
        main
        """
        metrics = self.read_quality_metrics()
        reads = self.read_corrected_reads()
        segs = self.read_segments()
        params = self.read_params()


        metrics = self.check_info_columns(metrics)

        plot_title = self.get_plot_title(self.sample_id, metrics)
        annotations = self.get_annotations(metrics)

        # If the check_mad returns false: filter it
        if not self.check_mad_score(self.sample_id, metrics):
            return

        # extract the data for the sample we're plotting
#         reads_samp = self.get_sample_data(reads, self.sample_id, norm=True)
#         segs_samp = self.get_sample_data(segs, self.sample_id)

        self.plot_corrected_reads(reads, self.sample_id, plot_title, annotations)
 
        self.plot_bias(reads, self.sample_id, plot_title, annotations)
 
        self.plot_segments(reads, segs, plot_title, annotations,
                           num_states=self.num_states)

        self.plot_params(reads, params, plot_title, annotations,
                         num_states=self.num_states)

        self.reads_pdf.close()
        self.bias_pdf.close()
        self.segs_pdf.close()
        self.params_pdf.close()

if __name__ == '__main__':
    args = parse_args()

    genhmm = GenHmmPlots(args.corrected_reads, args.segments, args.params, args.quality_metrics, args.sample_info, args.ref_genome,
                         args.reads_output, args.segs_output, args.bias_output, args.params_output, args.sample_id, annotation_cols=args.annotation_cols,
                         num_states = args.num_states, mad_threshold = args.mad_threshold, plot_title = args.plot_title)

    genhmm.main()
