'''
Created on Sep 8, 2015

@author: dgrewal
'''
import pandas
import matplotlib
matplotlib.use('Agg')
import seaborn as sns

from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import warnings

class PlotKernelDensity(object):
    '''
    merges files. no overlap queries, simple concatenation
    since columns are different, select header and insert values at proper
    indices. use N/A for missing.
    '''

    def __init__(self, infile, output,  sep, colname, plot_title, **kwargs):

        self.input = infile
        self.output = output
        self.sep = sep
        self.kernel_type = kwargs.get('kernel_type')
        self.column_name = colname
        self.plot_title = plot_title

        if not self.kernel_type:
            self.kernel_type = 'gau'

        bw_est = kwargs.get('bw_est')
        if not bw_est:
            bw_est = 'scott'
        
        if bw_est not in ["scott", "silverman"]:
            self.bw_est = float(bw_est)
        else:
            self.bw_est = bw_est

    def load(self, fname):
        '''
        load tsv file into a pandas data frame
        '''
        data = pandas.read_csv(fname,
                               sep=self.sep,
                               dtype={'chromosome': str, 'start': int}
                               )
        return data

    def plot_kernel_density(self, pdfout, mad_scores, exp_cond = None):
        """
        plots kernel density
        """
        sns.set_style('whitegrid')

        label = exp_cond if exp_cond else'all' 

        fig = sns.kdeplot(mad_scores,
                        bw=self.bw_est,
                        kernel=self.kernel_type,
                        label=label)

        fig.set(ylabel="Density", xlabel=self.column_name)

        fig = fig.get_figure()

        return fig

    def generate_plots(self):

        pdfout = PdfPages(self.output)
        data = self.load(self.input)

        #plot all data
        mad_scores = data[self.column_name]
        self.plot_kernel_density(pdfout, mad_scores)

        #get all experimental conditions
        exp_conds = set(data["experimental_condition"])
        # plot each exp cond
        for expcond in exp_conds:
            mad_scores = data[data["experimental_condition"] == expcond]\
                             [self.column_name]

            if mad_scores.isnull().all():
                warnings.warn("all mad states in condition %s are NaN" %expcond)
                continue

            fig = self.plot_kernel_density(pdfout, mad_scores, exp_cond=expcond)

        plt.suptitle(self.plot_title, fontsize=12)
        pdfout.savefig(fig)

        pdfout.close()


    def main(self):
        '''
        main function
        '''
        self.generate_plots()
