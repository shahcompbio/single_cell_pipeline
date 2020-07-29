'''
Created on Oct 24, 2017

@author: dgrewal
'''
import logging

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.cluster.hierarchy as hc
import scipy.spatial.distance as dist
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib.colors import rgb2hex


class ClusterMap(object):

    def __init__(self, data, colordata, max_cn, chromosomes=None,
                 scale_by_cells=False, distance_matrix=None):
        """
        :param data pandas dataframe with bins as columns and samples as rows
        :param colordata: dict with samples and their corresponding type
                        used for adding a colorbar
        """

        if chromosomes:
            self.chromosomes = chromosomes
        else:
            self.chromosomes = [str(v) for v in range(1, 23)] + ['X', 'Y']

        self.scale_by_cells = scale_by_cells

        self.max_cn = max_cn + 1

        self.colordata = colordata
        self.rows = data.index

        self.bins = data.columns.values

        self.data = data.as_matrix()

        # set max for data
        self.data = np.clip(self.data, 0, self.max_cn)

        self.distance_matrix = distance_matrix

        self.generate_plot()

    def get_chr_idxs(self, bins):
        """
        :param bins: sorted bins used for the plot
        :return chr_idxs: list with the index where chromosome changes
        returns the index where the chromosome changes
        used for marking chr boundaries on the plot
        """
        # chr 1 starts at beginning
        chr_idxs = [0]

        chrom = bins[0][0]
        for i, bin_v in enumerate(bins):
            if bin_v[0] != chrom:
                chr_idxs.append(i)
                chrom = bin_v[0]

        return chr_idxs

    def get_cmap_colorbar(self):
        """generates listed colormap for colordata
        used to plot a color bar. using seaborn to keep colors same as old code
        :returns listedcolormap
        :returns list of string desc for each color
        """
        ccs = list(set([self.colordata[row] for row in self.rows]))

        cmap = sns.color_palette("RdBu_d", len(ccs))

        return ListedColormap(cmap), ccs

    def generate_colormap_heatmap(self, maxval, vmin, vmax):
        """generating a custom heatmap 2:gray 0: blue 2+: reds
        :param maxval highest value in the data
        :returns listedcolormap
        """

        color_reference = {0: '#3498DB', 1: '#85C1E9', 2: '#D3D3D3'}

        low_max = 3 + int((maxval - 3) / 2) + 1
        hi_max = maxval
        low_states = np.arange(3, low_max)
        hi_states = np.arange(low_max, hi_max)

        low_cmap = matplotlib.cm.get_cmap('OrRd', low_max + 1)
        hi_cmap = matplotlib.cm.get_cmap('RdPu', hi_max + 1)

        for cn_level in low_states:
            rgb = low_cmap(int(cn_level))[:3]
            color_reference[cn_level] = rgb2hex(rgb)

        for cn_level in hi_states:
            rgb = hi_cmap(int(cn_level))[:3]
            color_reference[cn_level] = rgb2hex(rgb)

        color_reference[maxval] = "#000000"

        colors = [color_reference[cnlevel]
                  for cnlevel in np.arange(vmin, vmax + 1)]

        cmap = ListedColormap(colors)

        return cmap

    def plot_dendrogram(self, fig, mat, placement):
        """plots a dendrogram
        :param fig: matplotlib figure
        :param mat: data matrix (np.array)
        :param placement: list with [x,y,w,h] values for positioning plot
        :returns linkage matrix
        """

        # placement of dendrogram on the left of the heatmap
        ax1 = fig.add_axes(placement, frame_on=True)

        if self.distance_matrix is None:
            mat_no_nan = np.copy(mat)
            mat_no_nan[np.isnan(mat_no_nan)] = -1
            linkage = hc.linkage(dist.pdist(mat_no_nan, "cityblock"),
                                 method='ward')
        else:
            linkage = hc.linkage(self.distance_matrix, method='average')

        hc.dendrogram(linkage, orientation='left')
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.set_facecolor("white")

        return linkage

    def plot_row_colorbar(self, fig, linkage, placement):
        """adds colorbar next to the dendrogram
        :param fig: matplotlib figure
        :param linkage matrix
        :param placement: list with [x,y,w,h] values for positioning plot
        """
        axr = fig.add_axes(placement)

        order = hc.leaves_list(linkage)
        order = [self.rows[i] for i in order]

        cmap, ccs = self.get_cmap_colorbar()

        colors = np.array([ccs.index(self.colordata[v]) for v in order])
        colors.shape = (len(colors), 1)

        axr.matshow(colors, aspect='auto', origin='lower', cmap=cmap)
        axr.set_xticks([])
        axr.set_yticks([])

        # Plot color legend
        legend_plc = [0.77, 0.98, 0.18, 0.01]
        axcb = fig.add_axes(legend_plc, frame_on=False)

        self.plot_legend(axcb, cmap, ticklabels=ccs)

    def plot_heatmap(self, fig, mat, linkage, placement):
        """adds heatmap
        :param fig: matplotlib figure
        :param data matrix
        :param linkage matrix
        :param placement: list with [x,y,w,h] values for positioning plot
        """

        # sort matrix based on dendrogram order
        leaves = hc.leaves_list(linkage)
        mat = mat[leaves, :]

        if np.isnan(mat).all():
            logging.getLogger("single_cell.plot_heatmap").warn("skipping heatplot due to missing data")
            return

        # get all values we're going to plot later and generate a colormap for
        # them
        cmap = self.generate_colormap_heatmap(
            self.max_cn,
            np.nanmin(mat),
            np.nanmax(mat))

        # calculate appropriate margin to accomodate labels
        label_len = max([len(self.rows[leaves[i]])
                         for i in range(mat.shape[0])])
        right_margin = label_len * 0.005
        placement[-2] -= right_margin

        axm = fig.add_axes(placement)

        mat = np.ma.masked_where(np.isnan(mat), mat)

        axm.pcolormesh(
            mat,
            cmap=cmap,
            rasterized=True,
            vmin=np.nanmin(mat),
            vmax=np.nanmax(mat))

        axm.set_yticks([])

        for i in range(mat.shape[0]):
            axm.text(mat.shape[1] - 0.5, i, self.rows[leaves[i]],
                     fontsize=12)

        chr_idxs = self.get_chr_idxs(self.bins)
        axm.set_xticks(chr_idxs)
        axm.set_xticklabels(self.chromosomes)

        for val in chr_idxs:
            axm.plot(
                [val, val], [mat.shape[0], 0], ':', linewidth=0.5, color='black', )

        # Plot color legend
        legend_plc = [0.07, 0.98, 0.18, 0.01]
        axcb = fig.add_axes(legend_plc, frame_on=False)

        ticklabels = range(0, self.max_cn)
        ticklabels = list(map(str, ticklabels))
        ticklabels.append('{}+'.format(self.max_cn - 1))
        self.plot_legend(axcb, cmap, ticklabels=ticklabels)

    def plot_legend(self, axes, cmap, ticklabels=None):
        """adds legend
        :param axes: matplotlib figure axes
        :param cmap: colormap
        :param ticklabels: tick labels
        """

        bounds = range(cmap.N + 1)
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

        cbar = matplotlib.colorbar.ColorbarBase(axes, cmap=cmap, norm=norm,
                                                orientation='horizontal')
        cbar.set_ticks([v + 0.5 for v in bounds])

        if not ticklabels:
            ticklabels = [
                str(v).replace(str(self.max_cn), str(self.max_cn - 1) + "+") for v in bounds]

        if ticklabels:
            cbar.set_ticklabels(ticklabels)
        else:
            cbar.set_ticklabels(bounds)

    def generate_plot(self):
        """generates a figure with dendrogram, colorbar, heatmap and legends
        """

        figsize = (30, 30)
        if self.scale_by_cells:
            height = float(len(self.data)) / 7
            figsize = (30, height)

        fig = plt.figure(figsize=figsize)

        # fig's height and starting pos
        y = 0.1
        h = 0.85

        # dendrogram figure placement on the page
        dgram_plc = [0.05, y, 0.05, h]
        linkage = self.plot_dendrogram(fig, self.data, dgram_plc)

        # colorbar placement
        # x = dgram x + dgram w + margin
        # w=colorbar width
        cbar_plc = [0.101, y, 0.015, h]
        self.plot_row_colorbar(fig, linkage, cbar_plc)

        # heatmap placement
        # x = cbar x + cbar w + margin
        # w = width
        # right will get scaled based on labels inside plot_heatmap
        hmap_plc = [0.117, y, 0.9, h]
        self.plot_heatmap(fig, self.data, linkage, hmap_plc)

        return fig
