'''
Created on Nov 16, 2016

@author: dgrewal
'''

import logging
import os
from math import pi

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy.special import gamma


def t_dist_pdf(x, mu, lmbda, nu):
    p = (gamma(nu / 2 + 0.5) / gamma(nu / 2)) * \
        pow((lmbda / (pi * nu)), 0.5) * \
        pow(1 + (lmbda * pow((x - mu), 2)) / nu, -0.5 * nu - 0.5)
    return p


def read_fai(fai_filename):
    chrom_info = pd.read_csv(
        fai_filename, sep='\t', header=None,
        names=['chrom', 'length', 'V3', 'V4', 'V5'],
        usecols=['chrom', 'length'],
        dtype={'chrom': str, 'length': int})
    return chrom_info


def read_chromosome_lengths(ref_genome):
    chrom_info = read_fai(ref_genome + '.fai')
    chromosomes = list(chrom_info['chrom'])
    chrom_info = chrom_info[chrom_info['chrom'].isin(chromosomes)]
    chrom_length = chrom_info.set_index('chrom')['length']
    return chrom_length


def extract_chromosome_info(ref_genome):
    chromosome_info = pd.DataFrame()
    chromosome_info['lengths'] = read_chromosome_lengths(ref_genome)
    chromosome_info['end'] = np.cumsum(chromosome_info['lengths'])
    chromosome_info['start'] = chromosome_info['end'].shift(1)
    chromosome_info['start'][0] = 0
    chromosome_info['mid'] = (chromosome_info['start'] + chromosome_info['end']) / 2.
    return (chromosome_info)


def create_chromosome_plot_axes(ax, ref_genome):
    chromosome_info = extract_chromosome_info(ref_genome)
    ax.set_xlim((-0.5, chromosome_info['end'].max()))
    ax.set_xticks([0] + list(chromosome_info['end'].values))
    ax.set_xticklabels([])
    ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(chromosome_info['mid']))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.FixedFormatter(chromosome_info.index.values))
    ax.tick_params(which='minor', length=0)
    return (ax)


def compute_chromosome_coordinates(df, ref_genome):
    chromosome_info = extract_chromosome_info(ref_genome)
    df = df[df['chr'].isin(chromosome_info.index.values)]
    df.set_index('chr', inplace=True)
    df['chromosome_start'] = chromosome_info['start']
    df.reset_index(inplace=True)
    df['plot_coord'] = df['start'] + df['chromosome_start']
    df = df.drop('chromosome_start', 1)
    return (df)


def normalize_reads(df):
    if 'valid' not in df.columns.values:
        df['norm'] = float('nan')
    else:
        df['norm'] = df['reads'] / np.median(df['reads'])
        df['norm'] = df['norm'].where(df['valid'] == True)
    return (df)


def get_sample_id(out_file):
    sample_id = os.path.basename(out_file).split('.')[0]

    return sample_id


def add_legend(ax, color_reference, num_columns, type='rectangle', location='upper center'):
    object_list = []
    labels_list = []

    labels = sorted(color_reference.keys())
    colours = [color_reference[lbl] for lbl in color_reference]

    # legend for every fifth color if number is too high to make it fit in figure
    tick_size = 1
    if len(colours) > 25:
        tick_size = 5

    for i, col in enumerate(colours):
        if not i % tick_size == 0:
            continue

        labels_list.append(labels[i])

        if type == 'rectangle':
            object_list.append(Rectangle((0, 0), 1, 1, facecolor=col, edgecolor='none'))

        elif type == 'circle':
            object_list.append(Line2D(range(1), range(1), color=col, marker='o', markersize=15, linewidth=0))

        else:
            logging.getLogger("single_cell.plot_hmmcopy").warn(
                'Legend type must be one of: rectangle, circle.'
            )

    return ax.legend(tuple(object_list), tuple(labels_list), loc=location, ncol=num_columns)


def add_open_grid_lines(ax):
    for x_major in ax.xaxis.get_majorticklocs()[1:-1]:
        ax.axvline(x=x_major, linestyle=':', linewidth=1, color='.8', zorder=0)

    for y_major in ax.yaxis.get_majorticklocs()[1:-1]:
        ax.axhline(y=y_major, linestyle=':', linewidth=1, color='.8', zorder=0)


def get_segment_start_end(segments, remove_y=False):
    segment_diff = segments['plot_coord'] + (segments['end'] - segments['start'])
    segment_end = segments['plot_coord'][1:] - 1
    segment_end = segment_end.append(segment_diff.tail(1))
    segment_end.reset_index(inplace=True, drop=True)
    segments['plot_coord_end'] = segment_end

    if remove_y:
        segments = segments[segments['chr'] != 'Y']

    x = []
    y = []
    for x_start, x_end, y_med in zip(segments['plot_coord'], segments['plot_coord_end'], segments['median']):
        x.append(x_start)
        x.append(x_end)
        x.append(None)
        y.append(y_med)
        y.append(y_med)
        y.append(None)

    return (x, y)
