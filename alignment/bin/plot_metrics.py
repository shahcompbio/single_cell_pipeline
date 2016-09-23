'''
Plot sequencing metrics based on metric table and SampleSheet files.
'''
from __future__ import division

import argparse
import os
import matplotlib
matplotlib.use('Agg') # required for running on the cluster
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
matplotlib.rcParams['pdf.fonttype'] = 42

#=======================================================================================================================
# Read Command Line Input
#=======================================================================================================================
parser = argparse.ArgumentParser()

parser.add_argument('sample_sheet',
                    help='''Path to NextSeq sample sheet.''')

parser.add_argument('metric_table',
                    help='''Path to metric table for the run.''')

parser.add_argument('out_file',
                    help='''Path to output file where .pdf plots will be written.''')

args = parser.parse_args()

#=======================================================================================================================
# Functions
#=======================================================================================================================
def read_sample_sheet_data_table(file):
    with open(file) as f:
        lines = [x.strip('\n').strip(',') for x in f.readlines()]
    
    header = lines[lines.index('[Data]')+1].lower().split(',')
    
    data_lines = lines[lines.index('[Data]')+2:]
    
    df = pd.DataFrame(columns=header)
    
    for line in data_lines:
        sample_dict = dict(zip(header, line.split(',')))
        df = df.append(sample_dict, ignore_index=True)
    
    df = parse_sample_description(df)
    
    return df

def parse_sample_description(sample_sheet):
    '''
    Does not assume the description acronyms are in the same order in every sheet, 
    but does assume they are in the same order for every sample within a sheet.
    '''
    split = sample_sheet['description'].str.split(';')
    
    cell_call_column = np.where(['CC=' in x for x in split[0]])[0][0]
    exp_condition_column = np.where(['EC=' in x for x in split[0]])[0][0]
    
    sample_sheet['cell_call'] = [x[cell_call_column].replace('CC=', '') for x in split]
    
    sample_sheet['experimental_condition'] = [x[exp_condition_column].replace('EC=', '') for x in split]
    
    return sample_sheet

def add_legend(ax, labels, colours, num_columns, type='rectangle', location='upper center'):
    object_list = []
    for col in colours:
        if type=='rectangle':
            object_list.append(Rectangle((0, 0), 1, 1, facecolor=col, edgecolor='none'))
        
        elif type=='circle':
            object_list.append(Line2D(range(1), range(1), color=col, marker='o', markersize=15, linewidth=0))
        
        else:
            warnings.warn('Legend type must be one of: rectangle, circle.')
    
    return ax.legend(tuple(object_list), tuple(labels), loc=location, ncol=num_columns)

def add_barplot_labels(ax, labels, text_spacing, font_size):
    patch_width = list(set([x.get_width() for x in ax.patches]))[0]
    
    for plot_patch, label in zip(ax.patches, labels):
        patch_height = plot_patch.get_height()
        
        if np.isnan(patch_height):
            patch_height = 0
        
        ax.text(
                plot_patch.get_x() + patch_width / 2, 
                patch_height + text_spacing, 
                label,
                ha='center', 
                va='bottom',
                size=font_size
                )
    
    return ax

def plot_metric_fraction_total(df, metric, metric_label, pdf, from_top=False):
    sns.set(context='talk', 
            style='ticks', 
            font='Helvetica',
            rc={'axes.titlesize': 12,
                'axes.labelsize': 12, 
                'xtick.labelsize': 12,
                'ytick.labelsize': 12,
                'legend.fontsize': 12})
    
    fig = plt.figure(figsize=(len(df['sample_id'])/4, 5))
    
    ax = fig.gca()
    
    total = df['total_reads']/1000000
    fraction = df[metric]/1000000
    
    col_total = '#cfcfcf'
    col_fraction = '#595959'
    
    if from_top:
        fraction = total - fraction
        sns.barplot(df['sample_id'], total, color=col_fraction, ax=ax)
        sns.barplot(df['sample_id'], fraction, color=col_total, ax=ax)
    else:
        sns.barplot(df['sample_id'], total, color=col_total, ax=ax)
        sns.barplot(df['sample_id'], fraction, color=col_fraction, ax=ax)
    
    column_labels = [str(x) for x in df['cell_call']]
    
    ax = add_barplot_labels(ax, column_labels, 0.05, 12)
    
    ax.set_xlabel('Sample')
    ax.set_ylabel('Number of reads (millions)')
    sns.despine(offset=10, trim=True)
    
    sample_condition = [' (' + str(x) + ')' for x in df['experimental_condition']]
    sample_labels = [x + y for x, y in zip(df['sample_id'], sample_condition)]
    
    ax.set_xticklabels(sample_labels)
    
    plt.xticks(rotation=90)
    
    add_legend(ax, ['Total', metric_label], [col_total, col_fraction], 1, location='upper right')
    
    pdf.savefig(bbox_inches='tight', pad_inches=0.4)
    plt.close()

def plot_metric_fraction(df, numerator_metric, denominator_metric, ylab, pdf):
    sns.set(context='talk', 
            style='ticks', 
            font='Helvetica',
            rc={'axes.titlesize': 12,
                'axes.labelsize': 12, 
                'xtick.labelsize': 12,
                'ytick.labelsize': 12,
                'legend.fontsize': 12})
    
    fig = plt.figure(figsize=(len(df['sample_id'])/4, 5))
    
    ax = fig.gca()
    
    fraction = df[numerator_metric] / df[denominator_metric]
    
    col_fraction = '#595959'
    
    sns.barplot(df['sample_id'], fraction, color=col_fraction, ax=ax)
    
    column_labels = [str(x) for x in df['cell_call']]
    
    ax = add_barplot_labels(ax, column_labels, 0.015, 12)
    
    ax.set_xlabel('Sample')
    ax.set_ylabel(ylab)
    plt.ylim(0,1)
    sns.despine(offset=10, trim=True)
    
    sample_condition = [' (' + str(x) + ')' for x in df['experimental_condition']]
    sample_labels = [x + y for x, y in zip(df['sample_id'], sample_condition)]
    
    ax.set_xticklabels(sample_labels)
    
    plt.xticks(rotation=90)
    
    pdf.savefig(bbox_inches='tight', pad_inches=0.4)
    plt.close()

def plot_metric(df, metric, ylab, text_spacing, pdf):
    sns.set(context='talk', 
            style='ticks', 
            font='Helvetica',
            rc={'axes.titlesize': 12,
                'axes.labelsize': 12, 
                'xtick.labelsize': 12,
                'ytick.labelsize': 12,
                'legend.fontsize': 12})
    
    fig = plt.figure(figsize=(len(df['sample_id'])/4, 5))
    
    ax = fig.gca()
    
    col = '#595959'
    
    sns.barplot(df['sample_id'], df[metric], color=col, ax=ax)
    
    column_labels = [str(x) for x in df['cell_call']]
    
    ax = add_barplot_labels(ax, column_labels, text_spacing, 12)
    
    ax.set_xlabel('Sample')
    ax.set_ylabel(ylab)
    sns.despine(offset=10, trim=True)
    
    sample_condition = [' (' + str(x) + ')' for x in df['experimental_condition']]
    sample_labels = [x + y for x, y in zip(df['sample_id'], sample_condition)]
    
    ax.set_xticklabels(sample_labels)
    
    plt.xticks(rotation=90)
    
    pdf.savefig(bbox_inches='tight', pad_inches=0.4)
    plt.close()

def plot_metric_heatmap(df, metric, title, pdf, size=72):
    matrix = np.empty((size,size,))
    matrix[:] = np.nan
    
    well_labels = np.empty((size,size,))
    well_labels[:] = np.nan
    
    df['int_cell_call'] = [x.replace('NTC', '0').replace('C', '') for x in df['cell_call']]
    
    for i in range(len(df)):
        row_index = int(df.ix[i, 'sample_plate'].split('_')[0].replace('R', ''))
        col_index = int(df.ix[i, 'sample_plate'].split('_')[1].replace('C', ''))
        
        matrix_value = float(df.ix[i, metric])
        matrix[row_index-1, col_index-1] = matrix_value
        
        label_value = int(df.ix[i, 'int_cell_call'])
        well_labels[row_index-1, col_index-1] = label_value
    
    sns.set(context='talk', 
            style='darkgrid', 
            font='Helvetica',
            rc={'axes.titlesize': 9,
                'axes.labelsize': 6, 
                'xtick.labelsize': 6,
                'ytick.labelsize': 6,
                'legend.fontsize': 6})
    
    fig = plt.figure(figsize=(15,12))
    
    tick_labels = [x+1 for x in range(size)]
    
    # cheat to get well labels that are different from the colour value
    ax = sns.heatmap(well_labels, 
                     linewidths=0.6, 
                     square=True, 
                     annot=True, 
                     cmap=None, 
                     xticklabels=False, 
                     yticklabels=False, 
                     cbar=False, 
                     annot_kws={'size': 6}, 
                     alpha=0, 
                     fmt='.3g')
    
    for text in ax.texts:
        text.set_color('black')
    
    sns.heatmap(matrix, 
                xticklabels=tick_labels, 
                yticklabels=tick_labels, 
                linewidths=0.6, 
                square=True, 
                cbar=True)
    
    plt.title(title)
    
    pdf.savefig(bbox_inches='tight', pad_inches=0.4)
    
    plt.close()

def plot_metric_boxplot(df, metric, ylab, pdf):
    sns.set(context='talk', 
            style='ticks', 
            font='Helvetica',
            rc={'axes.titlesize': 12,
                'axes.labelsize': 12, 
                'xtick.labelsize': 12,
                'ytick.labelsize': 12,
                'legend.fontsize': 12})
    
    fig = plt.figure(figsize=(6, 6))
    
    ax = fig.gca()
    
    tableau_10_medium = ['#729ece', '#ff9e4a', '#67bf5c', '#ed665d', '#ad8bc9', '#a8786e', '#ed97ca', '#a2a2a2', '#cdcc5d', '#6dccda']
    
    num_conditions = len(df['experimental_condition'].unique())
    
    cols = tableau_10_medium[0:num_conditions]
    
    sns.boxplot(df['experimental_condition'], df[metric], palette=cols, ax=ax)
    
    num_libs = []
    
    for condition in df['experimental_condition'].unique():
        num_libs.append(len(df[df['experimental_condition']==condition]))
    
    condition_labels = [x + '\n(n=' + str(y) + ')' for x, y in zip(df['experimental_condition'].unique(), num_libs)]
    
    ax.set_xticklabels(condition_labels)
    
    ax.set_xlabel('Experimental condition')
    ax.set_ylabel(ylab)
    
    sns.despine(offset=10, trim=True)
    
    pdf.savefig(bbox_inches='tight', pad_inches=0.4)
    plt.close()

def plot_metric_factorplot(df, metric, ylab, pdf):
    df_melt = pd.melt(df, id_vars=['sample_id', 'experimental_condition', 'cell_call'], value_vars=[metric])
    
    experimental_conditions_ordered = sorted(df['experimental_condition'].unique())
    cell_calls_ordered = sorted(df['cell_call'].unique())
    
    num_cell_calls = len(cell_calls_ordered)
    tableau_10_medium = ['#729ece', '#ff9e4a', '#67bf5c', '#ed665d', '#ad8bc9', '#a8786e', '#ed97ca', '#a2a2a2', '#cdcc5d', '#6dccda']
    cols = tableau_10_medium[0:num_cell_calls]
    
    num_libs = []
    
    for condition in experimental_conditions_ordered:
        num_call = []
        
        for call in cell_calls_ordered:
            num_call.append(str(len(df[(df['experimental_condition']==condition) & (df['cell_call']==call)])))
        
        num_call = ', '.join(num_call)
        num_libs.append(num_call)
    
    condition_labels = [x + '\n(n=' + y + ')' for x, y in zip(experimental_conditions_ordered, num_libs)]    
    
    sns.set(context='talk', 
            style='ticks', 
            font='Helvetica',
            rc={'axes.titlesize': 12,
                'axes.labelsize': 12, 
                'xtick.labelsize': 12,
                'ytick.labelsize': 12,
                'legend.fontsize': 12})
    
    fig_height = 6
    fig_width = 8.5
    
    fig = plt.figure(figsize=(fig_width,fig_height))
    ax = fig.gca()
    
    g = sns.factorplot('experimental_condition', 
                        'value', 
                        'cell_call', 
                        df_melt, 
                        kind='box',
                        order=experimental_conditions_ordered,
                        hue_order=cell_calls_ordered,
                        palette=cols, 
                        legend=False, 
                        size=fig_height, 
                        aspect=fig_width/fig_height)
    
    plt.legend()
    
    g.set_axis_labels('Experimental condition', ylab)    
    g.set_xticklabels(condition_labels)
    
    sns.despine(offset=10, trim=True)
    
    pdf.savefig(bbox_inches='tight', pad_inches=0.4)
    plt.close()

#=======================================================================================================================
# Run script
#=======================================================================================================================
'''
args.sample_sheet = '/share/lustre/archive/single_cell_indexing/NextSeq/bcl/160115_NS500668_0072_AH5N5TAFXX/SampleSheet.csv'
args.metric_table = '/share/lustre/asteif/projects/single_cell_indexing/alignment/AH5N5TAFXX/lysis_conditions/metrics/summary/AH5N5TAFXX_lysis_conditions.metrics_table.csv'
args.out_file = '/share/lustre/asteif/projects/single_cell_indexing/test/nextseq_plots_test.pdf'
'''

def main():
    metrics = pd.read_csv(args.metric_table)
    
    samples = read_sample_sheet_data_table(args.sample_sheet)
    #samples['sample_id'] = [x.replace('NextSeq-', '') for x in samples['sample_id']]
    
    df = samples.merge(metrics, on='sample_id', how='left')
    
    with PdfPages(args.out_file) as pdf:
        plot_metric_fraction_total(df, 'total_mapped_reads', 'Mapped', pdf, from_top=False)
        plot_metric_fraction_total(df, 'total_duplicate_reads', 'Duplicates', pdf, from_top=True)
        plot_metric_fraction_total(df, 'total_properly_paired', 'Properly paired', pdf, from_top=False)
        
        plot_metric_fraction(df, 'total_mapped_reads', 'total_reads', 'Fraction mapped of total', pdf)
        plot_metric_fraction(df, 'total_duplicate_reads', 'total_mapped_reads', 'Fraction duplicates of mapped', pdf)
        plot_metric_fraction(df, 'total_properly_paired', 'total_mapped_reads', 'Fraction properly paired of mapped', pdf)
        
        plot_metric(df, 'coverage_depth', 'Coverage depth', 0.0015, pdf)
        plot_metric(df, 'coverage_breadth', 'Coverage breadth', 0.0015, pdf)
        plot_metric(df, 'mean_insert_size', 'Mean insert size', 0.05, pdf)
        plot_metric(df, 'median_insert_size', 'Median insert size', 0.05, pdf)
        
        plot_metric_heatmap(df, 'total_reads', 'Total reads', pdf)
        plot_metric_heatmap(df, 'percent_duplicate_reads', 'Percent duplicate reads', pdf)
        plot_metric_heatmap(df, 'coverage_depth', 'Coverage depth', pdf)
        plot_metric_heatmap(df, 'coverage_breadth', 'Coverage breadth', pdf)
        
        plot_metric_factorplot(df, 'total_reads', 'Total reads', pdf)
        plot_metric_factorplot(df, 'total_mapped_reads', 'Total mapped reads', pdf)
        plot_metric_factorplot(df, 'percent_duplicate_reads', 'Percent duplicate reads', pdf)
        plot_metric_factorplot(df, 'coverage_depth', 'Coverage depth', pdf)
        plot_metric_factorplot(df, 'coverage_breadth', 'Coverage breadth', pdf)
        plot_metric_factorplot(df, 'mean_insert_size', 'Mean insert size', pdf)
        plot_metric_factorplot(df, 'median_insert_size', 'Median insert size', pdf)
        plot_metric_factorplot(df, 'estimated_library_size', 'Picard estimated library size', pdf)

if __name__ == '__main__':
    main()
