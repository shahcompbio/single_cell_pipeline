'''
Created on Sep 8, 2015

@author: dgrewal
'''
import os
import sys
import math
import argparse
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import seaborn as sns
from collections import defaultdict
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import warnings


from heatmap import ClusterMap

sys.setrecursionlimit(2000)


class PlotPcolor(object):
    '''
    merges files. no overlap queries, simple concatenation
    since columns are different, select header and insert values at proper
    indices. use N/A for missing.
    '''

    def __init__(self, infile, metrics, output, **kwargs):
        self.input = infile
        self.metrics = metrics
        self.output = output

        if kwargs.get("chromosomes"):
            self.chromosomes = kwargs.get("chromosomes")
        else:
            self.chromosomes = [str(v) for v in range(1, 23)] + ['X', 'Y']

        self.sep = kwargs.get('separator')
        self.column_name = kwargs.get('column_name')
        self.cellcalls = kwargs.get('cellcalls')
        self.mad_thres = kwargs.get('mad_threshold')
        self.reads_thres = kwargs.get('numreads_threshold')
        self.median_hmmcopy_reads_per_bin_thres = kwargs.get(
            'median_hmmcopy_reads_per_bin_threshold')
        self.high_memory = kwargs.get('high_memory')
        self.plot_title = kwargs.get('plot_title')

        self.scale_by_cells = kwargs.get("scale_by_cells")

        self.color_by_col = kwargs.get('color_by_col')
        self.plot_by_col = kwargs.get('plot_by_col')

        self.mappability_threshold = kwargs.get('mappability_threshold')

        if not self.color_by_col:
            self.color_by_col = 'cell_call'

        if not self.plot_by_col:
            self.plot_by_col = 'all'

        if self.sep == 'comma':
            self.sep = ','
        elif self.sep == 'tab':
            self.sep = '\t'

        if not self.sep:
            self.sep = ','

        self.max_cn = kwargs.get("max_cn")

        if self.max_cn:
            self.max_cn = int(self.max_cn) + 1
        else:
            self.max_cn = 20

        self.multiplier = kwargs.get('multiplier')

    def build_label_indices(self, header):
        '''
        gets all the label cols from file and builds
        a dict with their indices as values
        '''
        if isinstance(header, str):
            header = header.strip().split(self.sep)

        lbl_idx = {val: i for i, val in enumerate(header)}

        return lbl_idx

    def sort_bins_hdf(self, bins):

        sortedbins = []

        bins = bins.groupby("chr")

        chromosomes = map(str, range(1, 23)) + ["X", "Y"]

        for chrom in chromosomes:
            chrom_bins = bins.get_group(chrom)
            chrom_bins = chrom_bins[["start", "end"]]

            chrom_bins = sorted(chrom_bins.values.tolist())

            chrom_bins = [(chrom, v[0], v[1]) for v in chrom_bins]

            sortedbins += chrom_bins

        return sortedbins

    def read_segs(self):

        extension = os.path.splitext(self.input)[-1]

        if extension in ['.hdf', '.h5']:
            return self.read_segs_hdf()
        else:
            return self.read_segs_csv()

    def read_segs_hdf(self):
        bins = None

        with pd.HDFStore(self.input, 'r') as reads_store:

            data = []

            for tableid in reads_store.keys():

                table = reads_store[tableid]

                if self.multiplier:
                    table_multiplier = int(tableid.split('/')[-1])
                    if not table_multiplier == self.multiplier:
                        continue


                if not bins:
                    bins = table[['chr', 'start', 'end']]

                    bins = self.sort_bins_hdf(bins)

                table["bin"] = list(zip(table.chr, table.start, table.end))

                table = table.pivot(
                    index='cell_id',
                    columns='bin',
                    values='state')

                table = table.sort_values(bins, axis=0)

                data.append(table)

        data = pd.concat(data)

        return data

    def read_segs_csv(self):
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

            if self.multiplier:
                if not int(line[idxs['multiplier']]) == self.multiplier:
                    continue

            chrom = line[idxs['chr']]
            start = int(line[idxs['start']])
            end = int(line[idxs['end']])

            seg = (chrom, start, end)

            if float(line[idxs["map"]]) <= self.mappability_threshold:
                val = float("nan")

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
        bins = self.sort_bins_csv(bins)

        data = self.conv_to_matrix(data, bins, samples)
        data = self.get_pandas_dataframe(data, bins)

        return data

    def sort_bins_csv(self, bins):
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

    def read_metrics_csv(self, cndata):
        """
        read the input file
        """

        samples = cndata.index

        data = {}
        numread_data = {}
        reads_per_bin_data = {}

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

            if self.multiplier:
                multiplier = int(line[idxs["multiplier"]])
                if not multiplier == self.multiplier:
                    continue

            val = line[idxs["mad_neutral_state"]]

            val = float('nan') if val == "NA" else float(val)

            ec = 'all' if sep_col == 'all' else line[idxs[sep_col]]

            cc = line[idxs[color_col]]

            numreads = int(line[idxs['total_mapped_reads']])

            reads_per_bin = line[idxs['median_hmmcopy_reads_per_bin']]

            reads_per_bin = 0 if reads_per_bin == "NA" else float(
                reads_per_bin)

            if self.cellcalls and cc not in self.cellcalls:
                continue

            numread_data[sample_id] = numreads
            data[sample_id] = val
            reads_per_bin_data[sample_id] = reads_per_bin

            colordata[sample_id] = cc
            sepdata[ec].append(sample_id)

        return data, sepdata, colordata, numread_data, reads_per_bin_data

    def read_metrics_hdf(self, cndata):

        samples = cndata.index

        data = []

        with pd.HDFStore(self.metrics, 'r') as metrics_store:
            for tableid in metrics_store:

                if self.multiplier:
                    table_multiplier = int(tableid.split('/')[-1])
                    if not table_multiplier == self.multiplier:
                        continue


                table = metrics_store[tableid]
                data.append(table)

        data = pd.concat(data)
        data = data.reset_index()

        plot_groups = {}
        color_groups = {}
        mad_data = {}
        numreads_data = {}
        reads_per_bin_data = {}

        for _, row in data.iterrows():
            cell = row["cell_id"]

            if cell not in samples:
                continue

            if self.plot_by_col == "all":
                plot_group_name = "all"
            else:
                plot_group_name = row[self.plot_by_col]

            color_group_name = row[self.color_by_col]

            mad_score = row["mad_neutral_state"]
            nreads = row["total_mapped_reads"]
            reads_per_bin = row["median_hmmcopy_reads_per_bin"]

            mad_data[cell] = mad_score
            numreads_data[cell] = nreads
            reads_per_bin_data[cell] = reads_per_bin

            if plot_group_name not in plot_groups:
                plot_groups[plot_group_name] = []
            plot_groups[plot_group_name].append(cell)

            color_groups[cell] = color_group_name

        return mad_data, plot_groups, color_groups, numreads_data, reads_per_bin_data

    def read_metrics(self, cndata):

        extension = os.path.splitext(self.metrics)[-1]

        if extension in ['.h5', '.hdf']:
            return self.read_metrics_hdf(cndata)
        else:
            return self.read_metrics_csv(cndata)

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

    def filter_data(
            self, data, ccdata, mad_scores, numreads_data, reads_per_bin):
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

        if self.median_hmmcopy_reads_per_bin_thres:
            samples = [samp for samp in samples
                       if reads_per_bin[samp] >= self.median_hmmcopy_reads_per_bin_thres]

        data = data.loc[samples]
        return data

    def plot_heatmap(self, data, ccdata, title, lims, pdfout):
        """
        generate heatmap, annotate and save

        """
        ClusterMap(
            data,
            ccdata,
            lims,
            self.max_cn,
            chromosomes=self.chromosomes,
            scale_by_cells=self.scale_by_cells)

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

            self.plot_heatmap(pltdata, colordata, title, lims, pdfout)

        if not self.output:
            return

        sns.set_style('whitegrid')
        sns.set(font_scale=1.5)

        pdfout = PdfPages(self.output)

        if not data.values.size:
            warnings.warn("no data to plot")
            return

        vmax = np.nanmax(data.values)
        vmin = np.nanmin(data.values)
        lims = (vmin, vmax)

        for sep, samples in sepdata.iteritems():

            num_samples = len(samples)

            samples = set(samples).intersection(set(data.index))

            if len(samples) < 2:
                continue

            if len(samples) > 1000 and not self.high_memory:
                warnings.warn('The output file will only plot 1000 cells per page,'
                              ' add --high_memory to override')

                samples = sorted(samples)
                # plot in groups of 1000
                sample_sets = [samples[x:x + 1000]
                               for x in range(0, len(samples), 1000)]
                for samples in sample_sets:

                    genplot(data, samples)
            else:
                genplot(data, samples)

        pdfout.close()

    def main(self):
        '''
        main function
        '''
        data = self.read_segs()

        if data.empty:
            warnings.warn("no data to plot")
            open(self.output, "w").close()
            return

        mad_scores, sepdata, colordata, numread_data, reads_per_bin = self.read_metrics(
            data)

        data = self.filter_data(
            data,
            colordata,
            mad_scores,
            numread_data,
            reads_per_bin)

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

    parser.add_argument('--plot_title',
                        help='''title for the plot''')

    parser.add_argument('--cellcalls',
                        nargs='*',
                        help='''list of the target cell types ''')

    parser.add_argument('--mad_threshold',
                        type=float,
                        default=None,
                        help='''all cells that have low MAD won't be plotted''')

    parser.add_argument('--numreads_threshold',
                        type=int,
                        default=None,
                        help='''all cells that have low MAD won't be plotted''')

    parser.add_argument('--median_hmmcopy_reads_per_bin_threshold',
                        type=int,
                        default=None,
                        help='''all cells that have low quality won't be plotted''')

    parser.add_argument('--mappability_threshold',
                        type=float,
                        default=0.9,
                        help='sets all cells with mappability under threshold to nan')

    parser.add_argument('--plot_by_col',
                        default='all',
                        help='''Column name to use for grouping the heatmaps''')

    parser.add_argument('--color_by_col',
                        default='cell_call',
                        help='''column name to use for coloring the side bar in heatmap''')

    parser.add_argument('--max_cn',
                        default=20,
                        help='''maximum copynumber to plot in heatmap''')

    parser.add_argument('--scale_by_cells',
                        default=False,
                        action="store_true",
                        help="scale the height of plot by number of cells")

    parser.add_argument('--high_memory',
                        action='store_true',
                        help='set this flag to override the default limit of 1000 cells'
                        ' per plot. The code will use more memory and the pdf file size'
                        ' will depend on number of cells')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    ARGS = parse_args()
    m = PlotPcolor(ARGS.input, ARGS.metrics, ARGS.output, column_name=ARGS.column_name,
                   cellcalls=ARGS.cellcalls, mad_threshold=ARGS.mad_threshold, numreads_threshold=ARGS.numreads_threshold,
                   median_hmmcopy_reads_per_bin_threshold=ARGS.median_hmmcopy_reads_per_bin_threshold, high_memory=ARGS.high_memory, plot_title=ARGS.plot_title,
                   color_by_col=ARGS.color_by_col, plot_by_col=ARGS.plot_by_col,
                   separator=ARGS.separator, max_cn=ARGS.max_cn, scale_by_cells=ARGS.scale_by_cells, mappability_threshold=ARGS.mappability_threshold)
    m.main()
