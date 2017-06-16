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
from collections import defaultdict

lowess = sm.nonparametric.lowess

matplotlib.rcParams['pdf.fonttype'] = 42

sns.set(context='talk', 
        style='ticks', 
        font='Helvetica',
        rc={'axes.titlesize': 12,
            'axes.labelsize': 15, 
            'xtick.labelsize': 15, 
            'ytick.labelsize': 15,
            'legend.fontsize': 15})


def parse_args():
    #=======================================================================================================================
    # Read Command Line Input
    #=======================================================================================================================
    parser = argparse.ArgumentParser()

    parser.add_argument('--corrected_reads',
                        required=True, 
                        help='''Path to HMMcopy corrected reads output .csv file.''')
    
    parser.add_argument('--segments',
                        required=True, 
                        help='''Path to HMMcopy segments output .csv file.''')
    
    parser.add_argument('--hmm_metrics',
                        required=True, 
                        help='''HMM quality metrics file for the run, with 'mad_neutral_state' column.''')

    parser.add_argument('--sample_info',
                        required=True, 
                        help='''Sample info .csv file.''')

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

    parser.add_argument('--samples',
                        nargs='*',
                        help='''title of the plots''')
    
    args = parser.parse_args()
    return args


class GenHmmPlots(object):
    """
    generate the reads, bias and segment plots
    """ 
    def __init__(self, args):
        self.args = args
        self.reads_pdf, self.segs_pdf, self.bias_pdf = self.get_pdf_handles()



    def load_data_pandas(self, infile):
        """
        
        """
        data = pd.read_csv(infile,
                           sep=',')

        # HACK!
        data['sample_id'] = data['sample_id'].apply(lambda a: a.replace('-', '_'))
        data = data.groupby('sample_id')
        
        return data


    def read_hmm_metrics(self):
        """
        
        """
        
        df = self.load_data_pandas(self.args.hmm_metrics)

        return df


    def read_sample_info(self):
        """

        """

        df = self.load_data_pandas(self.args.sample_info)
        
        return df


    def read_corrected_reads(self):
        """
        
        """
        
        df = self.load_data_pandas(self.args.corrected_reads)
        
        return df
    
    def read_segments(self):
        """
        
        """
        
        df = self.load_data_pandas(self.args.segments)
        
        return df

    def get_sample_ids(self, df, sample_info):
        """
        
        """
        samples = df.groups.keys()

        samdata = defaultdict(list)
        for samp in samples:
            ec = sample_info.get_group(samp)['experimental_condition'].iloc[0]
        
            samdata[ec].append(samp)
        
        sams = []
        
        for ec in sorted(samdata.keys()):
            sams.extend(sorted(samdata[ec]))
        
        return sams

    def get_pdf_handles(self):
        """
        
        """
        
        reads_pdf = PdfPages(self.args.reads_output)
        bias_pdf = PdfPages(self.args.bias_output)
        segs_pdf = PdfPages(self.args.segs_output)


        return reads_pdf, segs_pdf, bias_pdf

    def get_plot_title(self, sample_id, sample_info, metrics):
        """
        
        """
        if 'cell_call' in sample_info.get_group(sample_id):
            cellcall = sample_info.get_group(sample_id)['cell_call'].iloc[0]
        else:
            cellcall='NA'

        if 'experimental_condition' in sample_info.get_group(sample_id):
            cond = sample_info.get_group(sample_id)['experimental_condition'].iloc[0]
        else:
            cond='NA'

        mad = metrics.get_group(sample_id)['mad_neutral_state'].iloc[0]
        mad = str('%.3f' % mad)
        ni = metrics.get_group(sample_id)['MSRSI_non_integerness'].iloc[0]
        ni = str('%.3f' % ni)

        title_str = [sample_id , '(cell call', cellcall, ', condition',
                      cond, ', neutral MAD ', mad, ', MSRSI Non Integerness ',
                      ni, ')']

        title_str = ' '.join(title_str) + self.args.plot_title 
        return title_str

    def get_mad_score(self, sample_id, metrics):
        """
        """
        mad = metrics.get_group(sample_id)['mad_neutral_state'].iloc[0]
        return mad


    def gen_reads_plot(self, df, ax, typ = 'norm'):

        col = '#595959'

        ax = utl.create_chromosome_plot_axes(ax, self.args.ref_genome)

        #only attempt to plot if data is available
        if df is not None:
            plt.scatter(df['plot_coord'], df[typ], color=col, s=4)
            plt_lowess = lowess(df[typ], df['plot_coord'], frac=0.01, return_sorted=False)
            plt.plot(df['plot_coord'], plt_lowess, color='black', linewidth=1.2)

        if typ == 'norm':
            ax.set_ylabel('Normalized reads per bin')
        elif typ == 'cor.gc':
            ax.set_ylabel('GC corrected reads per bin')
        elif typ == 'cor.map':
            ax.set_ylabel('GC and mappability \n corrected reads per bin')
        ax.tick_params(axis='x', which='minor', pad=9.1)
        ax = utl.add_open_grid_lines(ax)


    def plot_corrected_reads(self, df, sample_id, title):
        """
        
        """
        fig = plt.figure(figsize=(15,12))
        
        plt.subplot(3, 1, 1)
        ax = fig.gca()
        self.gen_reads_plot(df, ax, typ='norm')
        ax.set_title(title)
        ax.set_xlabel('')
    
        plt.subplot(3, 1, 2)
        ax = fig.gca()
        self.gen_reads_plot(df, ax, typ='cor.gc')
        ax.set_xlabel('')
        
        plt.subplot(3, 1, 3)
        ax = fig.gca()
        self.gen_reads_plot(df, ax, typ='cor.map')
        
        sns.despine(offset=10, trim=True)
        plt.tight_layout()
        self.reads_pdf.savefig(fig, pad_inches=0.2)
        plt.close()

    def plot_bias(self, df, sample_id, title):
        """
        """
        df_ideal = df[df['ideal']==True]
        
        col = '#006ba4'
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(9, 9))
        
        ax1.set_ylabel('Read count')
        ax1.set_title('Uncorrected')
        if df is not None:
            ax1.scatter(df_ideal['gc'], df_ideal['reads'], edgecolors=col, facecolors='none', alpha=0.1)
        
        
        ax2.set_title('Uncorrected')
        if df is not None:
            ax2.scatter(df_ideal['map'], df_ideal['reads'], edgecolors=col, facecolors='none', alpha=0.1)
        
        ax3.set_xlabel('GC content')
        ax3.set_ylabel('Normalized read count')
        ax3.set_title('GC corrected')
        not_null = df_ideal['cor.gc'].notnull()
        if df is not None:
            ax3.scatter(df_ideal['gc'][not_null], df_ideal['cor.gc'][not_null], edgecolors=col, facecolors='none', alpha=0.1)
        
        ax4.set_xlabel('Mappability')
        ax4.set_title('GC and mappability corrected')
        not_null = df_ideal['cor.map'].notnull()
        if df is not None:
            ax4.scatter(df_ideal['map'][not_null], df_ideal['cor.map'][not_null], edgecolors=col, facecolors='none', alpha=0.1)

        fig.suptitle(title)
        sns.despine(offset=10, trim=True)
        
        plt.tight_layout(rect=(0,0,1,0.95))
        self.bias_pdf.savefig(fig, pad_inches=0.2)
        plt.close()

    def plot_segments(self, df, segments, plot_title, num_states=7, remove_y=False):
        if df is not None and remove_y:
            df = df[df['chr'] != 'Y']
    
        # standard: 15,4
        # SA501X3F xenograft heatmap: 20.4, 4
        fig = plt.figure(figsize=(15,4))
        ax = fig.gca()
    
        ax = utl.create_chromosome_plot_axes(ax, self.args.ref_genome)
        ax.set_title(plot_title)
        ax.set_ylabel('Copy number')
        
        segment_states = range(1, num_states+1)
        segment_labels = [str(x) for x in range(num_states)]
        segment_labels[-1] = segment_labels[-1] + ' or more'
        segment_colours = ['#006ba4', '#5f9ed1', '#ababab', '#ffbc79', '#ff800e', '#c85200', '#8c3900']
    
        #we pass None if we don't have data
        if df is not None and segments is not None:
            cols = df['state']
            cols = cols.replace(segment_states, segment_colours)
            cols = cols[~df['copy'].isnull()]
            
            plt.scatter(df['plot_coord'], df['integer_copy_scale'], facecolors=cols, edgecolors='none', s=4)
            
            x, y = utl.get_segment_start_end(segments, remove_y)
            plt.plot(x, y, color = 'black', linewidth=1)
        
        #ax.set_ylim((0, plt.ylim()[1]))
        ax.set_ylim((0, 14))
        
        sns.despine(offset=10, trim=True)
        ax.tick_params(axis='x', which='minor', pad=9.1)
        
        ax.legend = utl.add_legend(ax, segment_labels, segment_colours, num_states, type='rectangle', location='upper center')
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
        df = utl.compute_chromosome_coordinates(df, self.args.ref_genome)

        return df


    def check_mad_score(self, sample, metrics):
        """
        
        """
        mad = self.get_mad_score(sample, metrics)

        # if mad_threshold is set to nonzero.
        #zero is defaults and means mad_threshold is not set. so no filtering
        if self.args.mad_threshold:
            if math.isnan(mad):
                return False

            if mad > self.args.mad_threshold:
                return False
        return True



    def main(self):
        """
        main
        """
        sample_info = self.read_sample_info()
        metrics = self.read_hmm_metrics()
        reads = self.read_corrected_reads()
        segs = self.read_segments()

        if self.args.samples:
            samples = self.args.samples
        else:
            samples = self.get_sample_ids(reads, sample_info)

        for sample in samples:
            plot_title = self.get_plot_title(sample, metrics, sample_info)

            #If the check_mad returns false: filter it
            if not self.check_mad_score(sample, metrics):
                continue
            
            #extract the data for the sample we're plotting
            reads_samp = self.get_sample_data(reads, sample, norm=True)
            segs_samp = self.get_sample_data(segs, sample)

            self.plot_corrected_reads(reads_samp, sample, plot_title)
            
            self.plot_bias(reads_samp, sample, plot_title)
            
            self.plot_segments(reads_samp, segs_samp, plot_title,
                               num_states=self.args.num_states)

        self.reads_pdf.close()
        self.bias_pdf.close()
        self.segs_pdf.close()

if __name__ == '__main__':
    args = parse_args()
    
    genhmm = GenHmmPlots(args)

    genhmm.main()
