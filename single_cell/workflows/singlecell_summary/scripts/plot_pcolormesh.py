'''
Created on Sep 8, 2015

@author: dgrewal
'''
import sys
import math
import argparse
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc
from collections import defaultdict
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import warnings
import resource


from heatmap import ClusterMap

sys.setrecursionlimit(2000)


class PlotPcolor(object):
    '''
    merges files. no overlap queries, simple concatenation
    since columns are different, select header and insert values at proper
    indices. use N/A for missing.
    '''

    def __init__(self, infile, metrics, order_data, output, **kwargs):
        self.input = infile
        self.metrics = metrics
        self.order_data = order_data
        self.output = output

        self.chromosomes = [str(v) for v in range(1, 23)] + ['X', 'Y']

        self.sep = kwargs.get('sep')
        self.column_name = kwargs.get('colname')
        self.cellcalls = kwargs.get('cellcalls')
        self.mad_thres = kwargs.get('mad_thres')
        self.reads_thres = kwargs.get('reads_thres')
        self.high_memory = kwargs.get('high_memory')
        self.plot_title = kwargs.get('plot_title')


        self.color_by_col = kwargs.get('color_by_col')
        self.plot_by_col = kwargs.get('plot_by_col')

        if not self.color_by_col:
            self.color_by_col = 'cell_call'

        if not self.plot_by_col:
            self.plot_by_col = 'all'
            
        if not self.sep:
            self.sep = ','


    def build_label_indices(self, header):
        '''
        gets all the label cols from file and builds
        a dict with their indices as values
        '''
        if isinstance(header, str):
            header = header.strip().split(self.sep)

        lbl_idx = {val: i for i, val in enumerate(header)}

        return lbl_idx

    def read_segs(self):
        """
        read the input file
        """
        data = {}

        bins = {}

        freader = open(self.input)

        header = freader.readline()
        idxs = self.build_label_indices(header)

        for line in freader:
            line = line.strip().split(self.sep)

            sample_id = line[idxs['cell_id']]

            val = line[idxs[self.column_name]]

            val = float('nan') if val == "NA" else float(val)

            chrom = line[idxs['chr']]
            start = int(line[idxs['start']])
            end = int(line[idxs['end']])

            seg = (chrom, start, end)

            if chrom not in bins:
                bins[chrom] = set()
            bins[chrom].add((start, end))

            # just a sanity check, not required
            if sample_id in data and seg in data[sample_id]:
                raise Exception("repeated val")

            if sample_id not in data:
                data[sample_id] = {}

            data[sample_id][seg] = val

        samples = sorted(data.keys())
        return data, bins, samples

    def read_metrics(self, cndata):
        """
        read the input file
        """

        samples = cndata.index

        data = {}
        numread_data = {}

        sepdata = defaultdict(list)
        colordata = {}

        freader = open(self.metrics)

        header = freader.readline()
        idxs = self.build_label_indices(header)

        color_col = self.color_by_col
        sep_col = self.plot_by_col

        for line in freader:
            line = line.strip().split(self.sep)

            sample_id = line[idxs['cell_id']]

            # skip samples that are just na or inf
            if sample_id not in samples:
                continue

            val = line[idxs["mad_neutral_state"]]

            val = float('nan') if val == "NA" else float(val)

            ec = 'all' if sep_col=='all' else line[idxs[sep_col]]

            cc = line[idxs[color_col]]

            numreads = int(line[idxs['total_mapped_reads']])

            if self.cellcalls and cc not in self.cellcalls:
                continue

            numread_data[sample_id] = numreads
            data[sample_id] = val

            colordata[sample_id] = cc
            sepdata[ec].append(sample_id)

        return data, sepdata, colordata, numread_data

    def sort_bins(self, bins):
        """
        sort the bins based on genomic coords
        """
        assert set(self.chromosomes) == set(bins.keys())

        sort_bins = []
        for chrom in self.chromosomes:
            bin_vals = bins[chrom]

            bin_vals = sorted(bin_vals)

            bin_vals = [(chrom, bin_v[0], bin_v[1]) for bin_v in bin_vals]

            sort_bins += bin_vals

        return sort_bins

    def conv_to_matrix(self, data, bins, samples):
        """
        convert dict to numpy array
        """
        outdata = {}

        for sample in samples:
            cndata = [data[sample][bin_v] for bin_v in bins]

            # skip sample if all vals are nan or inf
            if np.isnan(cndata).all() or np.isinf(cndata).all():
                continue

            outdata[sample] = cndata

        return outdata

    def get_pandas_dataframe(self, data, bins):
        """
        convert array into dataframe
        provides an elegant way to annotate samples on the plot
        to remove nan rows and mask NA values
        only adds ~5s to runtime
        """
        df = pd.DataFrame(data)
        df = df.T
        df.columns = bins

        #remove cells that dont have any data
        df = df.dropna()
        return df

    @staticmethod
    def write_cluster_order(outfile, order):
        for i, samp in enumerate(order):
            outfile.write(','.join([samp, str(i)]) + '\n')

    @staticmethod
    def get_order(data):
        row_linkage = hc.linkage(sp.distance.pdist(data.values),
                                 method='average')
        order = hc.leaves_list(row_linkage)
 
        samps = data.index
        order = [samps[i] for i in order]
        return order

    def filter_data(self, data, ccdata, mad_scores, numreads_data):
        """
        remove samples that dont pass filtering thresholds
        """

        samples = data.index

        if self.cellcalls:
            samples = [samp for samp in samples
                       if ccdata[samp] in self.cellcalls]

        # remove samples over mad threshold
        if self.mad_thres:
            samples = [samp for samp in samples
                       if not math.isnan(mad_scores[samp])
                       and mad_scores[samp] <= self.mad_thres]

        # remove samples that have low num reads
        if self.reads_thres:
            samples = [samp for samp in samples
                       if numreads_data[samp] >= self.reads_thres]

        data = data.loc[samples]
        return data

    def get_cluster_order(self, data):
        """
        calculate distance matrix for clustering,
        get the ordering of the cells in the clustering
        dump order to file
        """
 
        if not self.order_data:
            return

        if all((self.cellcalls, self.mad_thres, self.reads_thres)):
            mad_scores, sepdata, colordata, numread_data = self.read_metrics(data)
            data = self.filter_data(data, colordata, mad_scores, numread_data)
        else:
            samples = list(set(data.index))
            sepdata = {'all':samples}


        outfile = open(self.order_data, 'w')
        
        outfile.write('cell_id,%s_heatmap_order\n' %self.plot_by_col)

        for _, samples in sepdata.iteritems():
            samples = set(samples).intersection(set(data.index))
            if len(samples) < 2:
                continue
            pltdata = data.loc[samples]
            order = self.get_order(pltdata)

            self.write_cluster_order(outfile, order)

        outfile.close()

    def plot_heatmap(self, data, ccdata, title, vmax, pdfout):
        """
        generate heatmap, annotate and save

        """
        ClusterMap(data, ccdata, vmax)

        plt.suptitle(title)

        plt.subplots_adjust(right=0.85)
 
        pdfout.savefig(pad_inches=0.2)

        plt.close("all")
        
    def plot_heatmap_by_sep(self, data, sepdata, colordata):
        """
        generate and save plot to output
        """
        def genplot(data, samples):
            pltdata = data.loc[samples]

            title = self.plot_title + \
                ' (%s) n=%s/%s' % (sep, len(samples), num_samples)

            self.plot_heatmap(pltdata, colordata, title, vmax, pdfout)

        
        if not self.output:
            return

        sns.set_style('whitegrid')
        sns.set(font_scale=1.5)

        pdfout = PdfPages(self.output)

        vmax = np.nanmax(data.values)


        for sep, samples in sepdata.iteritems():

            num_samples = len(samples)

            samples = set(samples).intersection(set(data.index))

            if len(samples) < 2:
                continue

            if len(samples) > 1000 and not self.high_memory:
                warnings.warn('The output file will only plot 1000 cells per page,'\
                              ' add --high_memory to override')

                samples = sorted(samples)
                #plot in groups of 1000
                sample_sets =  [samples[x:x+1000] for x in range(0, len(data), 1000)]
                for samples in sample_sets:

                    genplot(data, samples)
            else:
                genplot(data, samples)


        pdfout.close()

    def main(self):
        '''
        main function
        '''
        data, bins, samples = self.read_segs()

        bins = self.sort_bins(bins)
        data = self.conv_to_matrix(data, bins, samples)

        data = self.get_pandas_dataframe(data, bins)

        self.get_cluster_order(data)

        if self.output:
            mad_scores, sepdata, colordata, numread_data = self.read_metrics(data)
            data = self.filter_data(data, colordata, mad_scores, numread_data)

            self.plot_heatmap_by_sep(data, sepdata, colordata)



def parse_args():
    '''
    specify and parse args
    '''

    parser = argparse.ArgumentParser(description='''merge tsv/csv files''')

    parser.add_argument('--input',
                        required=True,
                        help='''corrected reads file from hmmcopy''')

    parser.add_argument('--metrics',
                        required=True,
                        help='''path to metrics file  ''')

    parser.add_argument('--column_name',
                        required=True,
                        help='column name of the value to be used'
                             ' for filling the values in heatmap')

    parser.add_argument('--separator',
                        required=True,
                        default="comma",
                        choices=("comma", "tab"),
                        help='''separator type, comma for csv, tab for tsv''')

    parser.add_argument('--output',
                        help='''path to output file''')

    parser.add_argument('--order_data',
                        help='''path to output clustering order''')

    parser.add_argument('--plot_title',
                        help='''title for the plot''')

    parser.add_argument('--cellcalls',
                        nargs='*',
                        help='''list of the target cell types ''')

    parser.add_argument('--mad_threshold',
                        type=float,
                        default=None,
                        dest='mad_thres',
                        help='''all cells that have low MAD won't be plotted''')

    parser.add_argument('--numreads_threshold',
                        type=int,
                        default=None,
                        dest='reads_thres',
                        help='''all cells that have low MAD won't be plotted''')

    parser.add_argument('--plot_by_col',
                        default='all',
                         help='''Column name to use for grouping the heatmaps''')

    parser.add_argument('--color_by_col',
                        default='cell_call',
                         help='''column name to use for coloring the side bar in heatmap''')

    parser.add_argument('--high_memory',
                        action='store_true',
                         help='set this flag to override the default limit of 1000 cells'\
                         ' per plot. The code will use more memory and the pdf file size'\
                         ' will depend on number of cells')


    args = parser.parse_args()

    return args


if __name__ == '__main__':
    ARGS = parse_args()
    m = PlotPcolor(ARGS.input, ARGS.metrics, ARGS.order_data, ARGS.output, colname=ARGS.column_name,
                    cellcalls=ARGS.cellcalls, mad_thres=ARGS.mad_thres, reads_thres=ARGS.reads_thres,
                    high_memory=ARGS.high_memory, plot_title=ARGS.plot_title, color_by_col=ARGS.color_by_col,
                    plot_by_col=ARGS.plot_by_col)
    m.main()
