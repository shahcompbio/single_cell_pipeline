'''
Created on Sep 8, 2015

@author: dgrewal
'''
import argparse
import pandas
import matplotlib
matplotlib.use('Agg')
import seaborn as sns

from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


class GenerateCNMatrix(object):
    '''
    merges files. no overlap queries, simple concatenation
    since columns are different, select header and insert values at proper
    indices. use N/A for missing.
    '''

    def __init__(self, args):
        self.args = args

        self.sep = ',' if self.args.separator == 'comma' else '\t'

        if self.args.bw_est not in ["scott", "silverman"]:
            self.bw_est = float(self.args.bw_est)
        else:
            self.bw_est = self.args.bw_est

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
                        kernel=self.args.kernel_type,
                        label=label)

        fig.set(ylabel="Density", xlabel=self.args.column_name)

        fig = fig.get_figure()

        return fig

    def generate_plots(self):

        pdfout = PdfPages(self.args.output)
        data = self.load(self.args.input)

        #plot all data
        mad_scores = data[self.args.column_name]
        self.plot_kernel_density(pdfout, mad_scores)

        #get all experimental conditions
        exp_conds = set(data["experimental_condition"])
        # plot each exp cond
        for expcond in exp_conds:
            mad_scores = data[data["experimental_condition"] == expcond]\
                             [self.args.column_name]
            fig = self.plot_kernel_density(pdfout, mad_scores, exp_cond=expcond)

        plt.suptitle(self.args.plot_title, fontsize=12)
        pdfout.savefig(fig)

        pdfout.close()


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
                        required=True,
                        default="comma",
                        choices=("comma", "tab"),
                        help='''separator type, comma for csv, tab for tsv''')

    parser.add_argument('--kernel_type',
                        default="gau",
                        choices=('gau', 'cos', 'biw', 'epa', 'tri', 'triw'),
                        help='''kernel type for the density estimation''')

    parser.add_argument('--bandwidth_estimator',
                        dest = "bw_est",
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
    m = GenerateCNMatrix(ARGS)
    m.main()
