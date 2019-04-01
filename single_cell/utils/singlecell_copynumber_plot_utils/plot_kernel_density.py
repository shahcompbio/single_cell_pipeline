'''
Created on Sep 8, 2015

@author: dgrewal
'''
import os
import pandas
import argparse
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import logging
import numpy as np
import statsmodels.nonparametric.api as smnp
import math

from single_cell.utils import csvutils

class PlotKernelDensity(object):
    '''
    merges files. no overlap queries, simple concatenation
    since columns are different, select header and insert values at proper
    indices. use N/A for missing.
    '''

    def __init__(self, infile, output, sep, colname, plot_title, **kwargs):

        self.input = infile
        self.output = output

        if sep == 'comma':
            self.sep = ','
        elif sep == 'tab':
            self.sep = '\t'
        else:
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

        self.tablename = kwargs.get("tablename")

    def load(self, fname):
        '''
        load tsv file into a pandas data frame
        '''
        extension = os.path.splitext(fname)[-1]

        if extension in [".h5", ".hdf5"]:

            with pandas.HDFStore(self.input, 'r') as metrics_store:
                data = metrics_store[self.tablename]

            data = data.reset_index()

        else:
            data = csvutils.read_csv_and_yaml(fname)

            # data['chromosome'] = data['chromosome'].astype(str)

        return data

    def get_ymax(self, data):
        if np.isnan(data).all():
            return 0

        kde = smnp.KDEUnivariate(data)
        kde.fit()

        maxval = np.nanmax(kde.density)
        if math.isnan(maxval):
            maxval = 0
        return maxval

    def plot_kernel_density(self, pdfout, mad_scores, exp_cond=None):
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

        return fig

    def generate_plots(self):

        pdfout = PdfPages(self.output)
        data = self.load(self.input)

        fig = None

        # plot all data
        mad_scores = data[self.column_name]

        if not np.isnan(mad_scores).all():
            self.plot_kernel_density(pdfout, mad_scores)
        ylim_max = self.get_ymax(mad_scores)

        # get all experimental conditions
        exp_conds = set(data["experimental_condition"])
        # plot each exp cond
        for expcond in exp_conds:
            mad_scores = data[data["experimental_condition"] == expcond]\
                             [self.column_name]

            if mad_scores.isnull().all():
                logging.getLogger("single_cell.plot_kernel_density").warn(
                    "all mad states in condition %s are NaN" %
                    expcond)
                continue

            if np.isnan(mad_scores).all():
                continue

            fig = self.plot_kernel_density(
                pdfout,
                mad_scores,
                exp_cond=expcond)

            ylim_max = max(ylim_max, self.get_ymax(mad_scores))

        if fig:
            fig.set_ylim((0, ylim_max))
            plt.tight_layout()
            plt.suptitle(self.plot_title, fontsize=12)
            pdfout.savefig(fig.get_figure())

        pdfout.close()
        plt.close()

    def main(self):
        '''
        main function
        '''
        self.generate_plots()


def parse_args():
    '''
    specify and parse args
    '''

    parser = argparse.ArgumentParser(
        description='''plot kernel density of MAD''')

    parser.add_argument('--input',
                        required=True,
                        help=''' path to the file with the qc data for all cells''')

    parser.add_argument('--separator',
                        default="comma",
                        choices=("comma", "tab"),
                        help='''separator type, comma for csv, tab for tsv''')

    parser.add_argument('--kernel_type',
                        default="gau",
                        choices=('gau', 'cos', 'biw', 'epa', 'tri', 'triw'),
                        help='''kernel type for the density estimation''')

    parser.add_argument('--bandwidth_estimator',
                        dest="bw_est",
                        default="scott",
                        help='''bandwith estimator type for picking the optimal bandwidth''')

    parser.add_argument('--column_name',
                        required=True,
                        help='''data from the specified column will be
                        extracted from input table and plotted''')

    parser.add_argument('--plot_title',
                        required=True,
                        help='''title for the plot''')

    parser.add_argument('--output',
                        required=True,
                        help='''path to output file''')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    ARGS = parse_args()
    kd = PlotKernelDensity(ARGS.input, ARGS.output, ARGS.separator,
                           ARGS.column_name, ARGS.plot_title,
                           kernel_type=ARGS.kernel_type,
                           bw_est=ARGS.bw_est)
    kd.main()
