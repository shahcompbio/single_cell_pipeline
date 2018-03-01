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
import matplotlib.gridspec as gridspec


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

    def load_data_pandas(self, infile, sample_id):
        """

        """
        data = pd.read_csv(infile,
                           sep=',')

        data = data[data['cell_id'] == sample_id]

        return data

    def load_data_pandas_lowmem(self, infile, sample_id, dtypes=None):
        chunksize = 10 ** 3
        data = pd.read_csv(infile, sep=',', chunksize=chunksize, dtype=dtypes)

        dfs= []
        for chunk in data:
            chunk = chunk[chunk["cell_id"] == sample_id]
            dfs.append(chunk)
        
        return pd.concat(dfs)


    def read_quality_metrics(self):
        """

        """

        df = self.load_data_pandas(self.metrics, self.sample_id)

        return df


    def read_params(self):
        """

        """

        df = self.load_data_pandas(self.params, self.sample_id)

        return df

    def read_corrected_reads(self):
        """

        """
        dtype = {"chr": str}

        df = self.load_data_pandas_lowmem(self.reads, self.sample_id, dtypes=dtype)

        df = utl.normalize_reads(df)
        df = utl.compute_chromosome_coordinates(df, self.ref_genome)

        return df

    def read_segments(self):
        """

        """
        dtype = {"chr": str}

        df = self.load_data_pandas_lowmem(self.segments, self.sample_id, dtypes=dtype)
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

    def add_annotations(self, fig, annotations, pos=(0.05,0.02), fontsize=None):
        annotations = ', '.join(annotations)

        fig.text(pos[0], pos[1], annotations, fontsize=fontsize)

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
            plt.scatter(df['plot_coord'], df[typ], color=col, s=4, rasterized=True)
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
        ax.tick_params(axis='x', which='minor', pad=9.1)
        ax = utl.add_open_grid_lines(ax)

    def plot_corrected_reads(self, df, sample_id, title, annotations):
        """

        """
        fig = plt.figure(figsize=(15, 12))

        plt.subplot(2, 1, 1)
        ax = fig.gca()
        self.gen_reads_plot(df, ax, typ='norm')
        ax.set_title(title)
        ax.set_xlabel('')

        plt.subplot(2, 1, 2)
        ax = fig.gca()
        self.gen_reads_plot(df, ax, typ='cor_gc')
        ax.set_xlabel('')

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

        fig = plt.figure(figsize=(12,15))

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

            
            mapp = {gc:pred for gc, pred in zip(df_ideal["gc"], df_ideal["modal_curve"])}
            x = sorted(df_ideal["gc"])
            y = [mapp[v] for v in x]
            plt.plot(x,y)


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

        self.add_annotations(fig, annotations, fontsize=10)
        fig.suptitle(title, fontsize=12)
        sns.despine(offset=10, trim=True)

        plt.tight_layout(rect=(0, 0.05, 1, 0.95))
        self.bias_pdf.savefig(fig, pad_inches=0.2)
        plt.close()

    def plot_segments(self, df, segments, params, plot_title,
                      annotations, num_states=7, remove_y=False):
        if df is not None and remove_y:
            df = df[df['chr'] != 'Y']

        ylim = num_states+1

        # standard: 15,4
        # SA501X3F xenograft heatmap: 20.4, 4
        fig = plt.figure(figsize=(18, 6))

        gs = gridspec.GridSpec(1, 2, width_ratios=[5, 1], wspace=0) 
        ax = plt.subplot(gs[0])
        
        ax = utl.create_chromosome_plot_axes(ax, self.ref_genome)
        ax.set_title(plot_title)
        ax.set_ylabel('Copy number')

        states = range(1, num_states+1)
        colors = ['#006ba4', '#5f9ed1', '#ababab', '#ffbc79', '#ff800e', '#c85200']
        #last color is same for all states > 6
        colors += (['#8c3900'] * (len(states) - len(colors)) )
        cmap = {i:v for i,v in zip(states, colors)}

        labels = [str(x) for x in range(num_states)]
        labels[-1] = labels[-1] + ' or more'

        # we pass None if we don't have data
        if df is not None and segments is not None:
            cols = df["state"].replace(cmap)
            cols = cols[~df['copy'].isnull()]
            
            plt.scatter(
                df['plot_coord'],
                df['integer_copy_scale'],
                facecolors=cols,
                edgecolors='none',
                s=4,
                rasterized=True)

            x, y = utl.get_segment_start_end(segments, remove_y)
            plt.plot(x, y, color='black', linewidth=1)

        #ax.set_ylim((0, plt.ylim()[1]))
        ax.set_ylim((0, ylim))

        sns.despine(offset=10, trim=True)
        ax.tick_params(axis='x', which='minor', pad=9.1)

        ax.legend = utl.add_legend(
            ax,
            labels,
            colors,
            num_states,
            type='rectangle',
            location='upper center')
        ax = utl.add_open_grid_lines(ax)

        self.add_annotations(fig, annotations)

        ax1 = plt.subplot(gs[1], sharey=ax)
        ax1.set_ylim((0, ylim))

        sns.despine(offset=10)

        ax1.yaxis.set_visible(False)
        ax1.spines['left'].set_visible(False)

        self.plot_dist(ax1, df, params, plot_title, annotations, cmap, num_states, vertical=True)

        plt.tight_layout(rect=(0, 0.05, 1, 1))
        self.segs_pdf.savefig(fig, pad_inches=0.2)
        plt.close()


    def plot_dist(self, fig, reads, params, plot_title,
                      annotations, cmap, num_states=7, vertical=False):


        #no data to plot
        if reads['integer_copy_scale'].isnull().all():
            return

        scale = (reads["integer_copy_scale"]/reads["copy"])
        scale = scale[~reads['copy'].isnull()].unique()
        #account for floating point errors
        assert np.nanvar(scale) < 0.0001
        scale = scale[0]

        for state in range(1, num_states+1):
            color = cmap[state]

            data = reads[reads["state"] == state]["integer_copy_scale"]


            if np.isnan(data).all():
                continue

            if not data.empty:
                sns.kdeplot(data, bw="scott",  kernel="epa",
                            shade=True, linewidth=0.5,
                            label=state, legend=False,
                            color = color, vertical=vertical)

            x = np.arange(0, np.nanmax(np.array(reads["copy"])), 0.01)
            mu = params[(params["parameter"]=="mus") & (params["state"] == state)]["final"].iloc[0]
            lmbda = params[(params["parameter"]=="lambdas") & (params["state"] == state)]["final"].iloc[0]
            nu = params[(params["parameter"]=="nus") & (params["state"] == state)]["final"].iloc[0]

            y = utl.t_dist_pdf(x, mu, lmbda, nu)
            x = x*scale

            if vertical:
                plt.plot(y,x, color, linewidth=0.5, linestyle='--')
                fig.set_xlim(0,4)
            else:
                plt.plot(x, y, color, linewidth=0.5, linestyle='--')
                fig.set_ylim(0,4)

        fig.set_xlabel("density")


    def plot_params(self, reads, params, plot_title,
                      annotations, num_states=7):

        plt.figure(figsize=(10, 20))
        
        fig,ax = plt.subplots()
        
        states = range(1, num_states+1)
        colors = ['#006ba4', '#5f9ed1', '#ababab', '#ffbc79', '#ff800e', '#c85200']
        #last color is same for all states > 6
        colors += (['#8c3900'] * (len(states) - len(colors)) )
        cmap = {i:v for i,v in zip(states, colors)}

        labels = [str(x) for x in range(num_states)]
        labels[-1] = labels[-1] + ' or more'

        self.plot_dist(ax, reads, params, plot_title, annotations, cmap, num_states)
        
        plt.suptitle(plot_title)

        self.add_legend(ax, labels, colors)

        plt.tight_layout(rect=(0, 0.05, 1, 0.95))
        self.add_annotations(ax, annotations, fontsize=9, pos=(-2, -0.6) )
        self.params_pdf.savefig(fig)
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

    def add_legend(self, fig, labels, colors):
        lgnd_patches = [Patch(color=c, label=k)
                        for k,c in zip(labels, colors)]
        plt.legend(handles=lgnd_patches,
                   bbox_to_anchor=(1, 1))


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

        self.plot_corrected_reads(reads, self.sample_id, plot_title, annotations)
  
        self.plot_bias(reads, self.sample_id, plot_title, annotations)
 
        self.plot_segments(reads, segs, params, plot_title, annotations,
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
