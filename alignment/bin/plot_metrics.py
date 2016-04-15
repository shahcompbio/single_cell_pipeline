'''
Plot sequencing metrics based on metric table and SampleSheet files.
'''
from __future__ import division

import argparse
import os
import matplotlib
#matplotlib.use('Agg') # required for running on the cluster
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

parser.add_argument('--sample_sheet',
                    help='''Path to NextSeq sample sheet.''')

parser.add_argument('--metric_table',
                    help='''Path to metric table for the run.''')

parser.add_argument('--out_file',
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
    
    return df

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



#=======================================================================================================================
# Run script
#=======================================================================================================================
'''
args.sample_sheet = '/share/lustre/archive/single_cell_indexing/NextSeq/bcl/160115_NS500668_0072_AH5N5TAFXX/SampleSheet.csv'
args.metric_table = '/share/lustre/asteif/projects/single_cell_indexing/alignment/AH5N5TAFXX/lysis_conditions/metrics/summary/AH5N5TAFXX_lysis_conditions.metrics_table.csv'
args.out_file = '/share/lustre/asteif/projects/single_cell_indexing/test/nextseq_plots_test.pdf'

Note: sample_well should be the location on the chip
'''


def plot_metric_heatmap(df, metric, size=72):

matrix = np.empty((size,size,))
matrix[:] = np.nan

for i in range(len(df)):
    row_index = int(df.ix[i, 'sample_well'].split('_')[0].replace('R', ''))
    col_index = int(df.ix[i, 'sample_well'].split('_')[1].replace('C', ''))
    value = float(df.ix[i, metric])
    matrix[row_index-1, col_index-1] = value

labels = np.empty((size,size,))
labels[:] = np.nan

for i in range(len(df)):
    row_index = int(df.ix[i, 'sample_well'].split('_')[0].replace('R', ''))
    col_index = int(df.ix[i, 'sample_well'].split('_')[1].replace('C', ''))
    value = str(df.ix[i, 'description'])
    labels[row_index-1, col_index-1] = value



# cheat to get labels that are different from the colour value!
ax = sns.heatmap(labels, linewidths=0.2, annot=True, cmap=None, xticklabels=False, yticklabels=False, cbar=False, annot_kws={'size': 8}, alpha=0)
for text in ax.texts:
    text.set_color('black')

sns.heatmap(matrix, linewidths=0.2, cbar=False)




def main():
metrics = pd.read_csv(args.metric_table)

samples = read_sample_sheet_data_table(args.sample_sheet)
#samples['sample_id'] = [x.replace('NextSeq-', '') for x in samples['sample_id']]

df = samples.merge(metrics, on='sample_id', how='left')

sns.set(context='talk', 
        style='darkgrid', 
        font='Helvetica',
        rc={'axes.titlesize': 8,
            'axes.labelsize': 8, 
            'xtick.labelsize': 8,
            'ytick.labelsize': 8,
            'legend.fontsize': 8})

if __name__ == '__main__':
    main()
