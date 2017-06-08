'''
Created on Sep 8, 2015

@author: dgrewal
'''
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import pyplot as plt
from matplotlib.colors import rgb2hex

import seaborn as sns
from matplotlib.colors import ListedColormap
from collections import defaultdict
import math


class PlotHeatmap(object):
    '''
    merges files. no overlap queries, simple concatenation
    since columns are different, select header and insert values at proper
    indices. use N/A for missing.
    '''

    def __init__(self, args):
        self.args = args

        self.sep = ',' if self.args.separator == 'comma' else '\t'

        self.chromosomes = [str(v) for v in range(1, 23)] + ['X', 'Y']

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

        freader = open(self.args.input)

        header = freader.readline()
        idxs = self.build_label_indices(header)

        for line in freader:
            line = line.strip().split(self.sep)

            sample_id = line[idxs['cell_id']]

            val = line[idxs[self.args.column_name]]

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
        ccdata = {}
        numread_data = {}

        ecdata = defaultdict(list)

        freader = open(self.args.metrics)

        header = freader.readline()
        idxs = self.build_label_indices(header)

        for line in freader:
            line = line.strip().split(self.sep)

            sample_id = line[idxs['cell_id']]

            #skip samples that are just na or inf
            if sample_id not in samples:
                continue

            val = line[idxs["mad_neutral_state"]]

            val = float('nan') if val == "NA" else float(val)

            ec = line[idxs["experimental_condition"]]

            if self.args.by_sample_type:
                ec = line[idxs["sample_type"]]

            cc = line[idxs["cell_call"]]

            numreads = int(line[idxs['total_mapped_reads']])
            
            
            if self.args.cellcalls and cc not in self.args.cellcalls:
                continue
            
            # just a sanity check, not required
            if sample_id in data:
                raise Exception("repeated val")

            # just a sanity check, not required
            if sample_id in ecdata:
                raise Exception("repeated val")

            numread_data[sample_id] = numreads
            data[sample_id] = val
            ccdata[sample_id] = cc
            ecdata[ec].append(sample_id)
        return data, ecdata, ccdata, numread_data

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

        return df

    def get_chr_idxs(self, bins):
        """
        returns the index where the chromosome changes
        used for marking chr boundaries on the plot
        """
        # chr 1 starts at beginning
        chr_idxs = [0]

        chrom = '1'
        for i, bin_v in enumerate(bins):
            if bin_v[0] != chrom:
                chr_idxs.append(i)
                chrom = bin_v[0]

        return chr_idxs

    def generate_colormap(self, maxval):
        """
        generating a custom heatmap 2:gray 0: blue 2+: reds
        """
        if self.args.column_name != 'integer_copy_number':
            return matplotlib.cm.coolwarm
        
        # all colors 2 and up are red with increasing intensity
        num_reds = maxval

        cmap = matplotlib.cm.get_cmap('Reds', num_reds)

        reds_hex = []
        for i in range(2, cmap.N):
            # will return rgba, we take only first 3 so we get rgb
            rgb = cmap(i)[:3]
            reds_hex.append(rgb2hex(rgb))

        cmap = ListedColormap(['#3498DB', '#85C1E9', '#D3D3D3'] + reds_hex)

        return cmap

    def get_colors(self, ccdata):
        """
        """
        if self.args.cellcalls:
            ccs = self.args.cellcalls
        else:
            ccs = list(set(ccdata.values()))
        colmap = sns.color_palette("RdBu_d", len(ccs))

        colmap = {cc:col for cc, col in zip(ccs, colmap)}

        return colmap

    def plot_heatmap(self, data, chr_idxs, mad_scores, ecdata, ccdata, numreads_data):
        """
        generate and save plot to output
        """
        sns.set_style('whitegrid')
        sns.set(font_scale=1.5)

        pdfout = PdfPages(self.args.output)
        cmap = self.generate_colormap(np.nanmax(data.values))
        vmax=np.nanmax(data.values)
        rowclr = self.get_colors(ccdata)

        for ec, samples in ecdata.iteritems():

            num_samples = len(samples)
            
            if self.args.cellcalls:
                samples = [samp for samp in samples
                           if ccdata[samp] in self.args.cellcalls]
            
            #remove samples over mad threshold
            if self.args.mad_thres:
                samples = [samp for samp in samples
                           if not math.isnan(mad_scores[samp])
                             and mad_scores[samp] <= self.args.mad_thres]
            
            #remove samples that have low num reads
            if self.args.reads_thres:
                samples = [samp for samp in samples
                           if numreads_data[samp] >= self.args.reads_thres]
            
           
            if len(samples) < 2:
                continue
 
            ecdata = data.loc[samples]
         
            mask = ecdata.isnull()

            colors = [rowclr[ccdata[samp]] for samp in samples]

            heatmap = sns.clustermap(ecdata, rasterized=True, mask=mask,
                                    figsize=(30, 50), cmap=cmap,
                                    vmin=0, vmax=vmax,
                                    col_cluster=False,
                                    row_colors = colors)
 
            heatmap.ax_heatmap.set(xticks=chr_idxs)
            heatmap.ax_heatmap.set(xticklabels=self.chromosomes)

            title=self.args.plot_title + ' (%s) n=%s/%s'%(ec, len(samples), num_samples)
            heatmap.ax_heatmap.set(title=self.args.plot_title + ' (%s)'%ec)

            plt.setp(heatmap.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

            # Plot the legend
            lgnd_patches = [matplotlib.patches.Patch(color=rowclr[k], label=k) for k in list(set(ccdata.values()))]
            heatmap.ax_heatmap.legend(handles=lgnd_patches, bbox_to_anchor=(1, 1.2))

            pdfout.savefig(pad_inches=0.2)
        pdfout.close()

    def main(self):
        '''
        main function
        '''
        data, bins, samples = self.read_segs()
        
        bins = self.sort_bins(bins)
        data = self.conv_to_matrix(data, bins, samples)

        data = self.get_pandas_dataframe(data, bins)
        chr_idxs = self.get_chr_idxs(bins)

        
        mad_scores, ecdata, ccdata, numread_data = self.read_metrics(data)

        self.plot_heatmap(data, chr_idxs, mad_scores, ecdata, ccdata, numread_data)


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
                        help='''column name of the value to be used for filling the values in heatmap''')

    parser.add_argument('--separator',
                        required=True,
                        default="comma",
                        choices=("comma", "tab"),
                        help='''separator type, comma for csv, tab for tsv''')

    parser.add_argument('--output',
                        required=True,
                        help='''path to output file''')

    parser.add_argument('--plot_title',
                        help='''title for the plot''')

    parser.add_argument('--cellcalls',
                        nargs='*',
                        help='''list of the target cell types ''')

    parser.add_argument('--by_sample_type',
                        default=False,
                        action='store_true',
                        help='''split plots by sample type instead of experimental condition ''')

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

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    ARGS = parse_args()
    m = PlotHeatmap(ARGS)
    m.main()

