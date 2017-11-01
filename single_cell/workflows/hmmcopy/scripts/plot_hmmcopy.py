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

lowess = sm.nonparametric.lowess

matplotlib.rcParams['pdf.fonttype'] = 42

sns.set(context='talk',
        style='ticks',
        font='Helvetica',
        rc={'axes.titlesize': 10,
            'axes.labelsize': 12,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'legend.fontsize': 12})


def parse_args():
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

    parser.add_argument('--plot_title',
                        help='''title of the plots''')

    parser.add_argument('--sample_id',
                        help='''title of the plots''')

    args = parser.parse_args()
    return args


class GenHmmPlots(object):
    """
    generate the reads, bias and segment plots
    """

    def __init__(self, reads, segments, metrics, sample_info,
                 ref_genome, reads_out, segs_out, bias_out, sample_id, **kwargs):

        self.reads = reads
        self.segments = segments
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

        self.reads_pdf, self.segs_pdf, self.bias_pdf = self.get_pdf_handles(
            reads_out, bias_out, segs_out)

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

    def get_pdf_handles(self, reads_out, bias_out, segs_out):
        """

        """

        reads_pdf = PdfPages(reads_out)
        bias_pdf = PdfPages(bias_out)
        segs_pdf = PdfPages(segs_out)

        return reads_pdf, segs_pdf, bias_pdf

    def get_plot_title(self, sample_id, metrics):
        """

        """
        if 'cell_call' in metrics:
            cellcall = metrics['cell_call'].iloc[0]
        else:
            cellcall = 'NA'

        if 'experimental_condition' in metrics:
            cond = metrics['experimental_condition'].iloc[0]
        else:
            cond = 'NA'

        if 'sample_type' in metrics:
            st = metrics['sample_type'].iloc[0]
            st = str(st)
        else:
            st = 'NA'

        mad = metrics['mad_neutral_state'].iloc[0]
        mad = str('%.3f' % mad)
        ni = metrics['MSRSI_non_integerness'].iloc[0]
        ni = str('%.3f' % ni)

        title_str = [sample_id, '(cell call', cellcall, ', condition',
                     cond, ', sample_type', st, ', neutral MAD ', mad, ', MSRSI NonInt. ',
                     ni, ')']

        title_str = ' '.join(title_str) + self.plot_title
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

    def plot_corrected_reads(self, df, sample_id, title):
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
        plt.tight_layout()
        self.reads_pdf.savefig(fig, pad_inches=0.2)
        plt.close()

    def plot_bias(self, df, sample_id, title):
        """
        """
        df_ideal = df[df['ideal'] == True]

        col = '#006ba4'

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
            2, 2, sharex='col', sharey='row', figsize=(9, 9))

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

        fig.suptitle(title, fontsize=10)
        sns.despine(offset=10, trim=True)

        plt.tight_layout(rect=(0, 0, 1, 0.95))
        self.bias_pdf.savefig(fig, pad_inches=0.2)
        plt.close()

    def plot_segments(
            self, df, segments, plot_title, num_states=7, remove_y=False):
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

        plt.tight_layout()
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

    def main(self):
        """
        main
        """
        metrics = self.read_quality_metrics()
        reads = self.read_corrected_reads()
        segs = self.read_segments()

        metrics = self.check_info_columns(metrics)

        plot_title = self.get_plot_title(self.sample_id, metrics)

        # If the check_mad returns false: filter it
        if not self.check_mad_score(self.sample_id, metrics):
            return

        # extract the data for the sample we're plotting
#         reads_samp = self.get_sample_data(reads, self.sample_id, norm=True)
#         segs_samp = self.get_sample_data(segs, self.sample_id)

        self.plot_corrected_reads(reads, self.sample_id, plot_title)

        self.plot_bias(reads, self.sample_id, plot_title)

        self.plot_segments(reads, segs, plot_title,
                           num_states=self.num_states)

        self.reads_pdf.close()
        self.bias_pdf.close()
        self.segs_pdf.close()

if __name__ == '__main__':
    args = parse_args()

    genhmm = GenHmmPlots(args.corrected_reads, args.segments, args.quality_metrics, args.sample_info, args.ref_genome,
                         args.reads_output, args.segs_output, args.bias_output, args.sample_id)

    genhmm.main()
