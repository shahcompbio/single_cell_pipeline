'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pandas as pd
from scripts import SummaryMetrics
from scripts import PlotKernelDensity
from scripts import PlotHeatmap
from scripts import PlotMetrics
from scripts import PlotPcolor
from scripts import classifier
from scripts import GenHmmPlots
from single_cell.utils import csvutils
from single_cell.utils import pdfutils
from single_cell.utils import helpers


def classify(metrics, output):

    scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
    model = os.path.join(scripts_directory, 'model.npz')
    classifier.main(metrics, model, output)

def annotate_metrics(infile, sample_info, outfile):
    infile = pd.read_csv(infile)

    cols = ["cell_call", "experimental_condition", "sample_type", "i5_barcode", "i7_barcode", "sample_plate", "sample_well"]

    for col in cols:
        infile[col] = "NA"

    cells = infile["cell_id"]
    
    for col in cols:
        coldata = [sample_info[cell][col] for cell in cells]
        infile[col] = coldata

    infile.to_csv(outfile, index=False, na_rep="NA")

def merge_tables(infile, output, sep,
                 merge_type, key_cols, nan_val):
    csvutils.merge_csv(infile, output, merge_type, key_cols)


def concatenate_csv(in_filenames, out_filename, nan_val='NA'):
    csvutils.concatenate_csv(in_filenames, out_filename)

def get_summary_metrics(infile, output):
    summ = SummaryMetrics(infile, output)
    summ.main()


def plot_kernel_density(infile, output, sep, colname, plot_title):
    plot = PlotKernelDensity(infile, output, sep, colname, plot_title)
    plot.main()

def plot_metrics(metrics, output, plot_title, gcbias_matrix, gc_content):
    plot = PlotMetrics(metrics, output, plot_title, gcbias_matrix, gc_content)
    plot.main()


def plot_heatmap(infile, metrics, order_data, output, plot_title=None,
                 column_name=None, plot_by_col=None, numreads_threshold=None,
                 mad_threshold=None, chromosomes=None, max_cn=None):

    plot = PlotHeatmap(infile, metrics, order_data, output, plot_title=plot_title,
                       column_name=column_name, plot_by_col=plot_by_col,
                       numreads_threshold=numreads_threshold,
                       mad_threshold=mad_threshold, chromosomes=chromosomes,
                       max_cn = max_cn)
    plot.main()

def plot_pcolor(infile, metrics, order_data, output, plot_title=None,
                 column_name=None, plot_by_col=None, numreads_threshold=None,
                 mad_threshold=None, chromosomes=None, max_cn=None, quality_threshold=None):

    plot = PlotPcolor(infile, metrics, order_data, output, plot_title=plot_title,
                       column_name=column_name, plot_by_col=plot_by_col,
                       numreads_threshold=numreads_threshold,
                       mad_threshold=mad_threshold, chromosomes=chromosomes,
                       max_cn = max_cn, quality_threshold = quality_threshold)
    plot.main()


def plot_hmmcopy(reads, segments, params, metrics, sample_info, ref_genome, reads_out, segs_out,
                 bias_out, params_out, sample_id, num_states=7, plot_title=None,
                 mad_threshold=None, annotation_cols=None):
    plot = GenHmmPlots(reads, segments, params, metrics, sample_info, ref_genome, reads_out, segs_out,
                       bias_out, params_out, sample_id, num_states=num_states, plot_title=plot_title,
                       mad_threshold=mad_threshold, annotation_cols=annotation_cols)
    plot.main()


def get_mad_score(samp, metrics):
    mad = metrics[metrics['cell_id'] == samp]['mad_neutral_state'].iloc[0]
    return mad

def get_num_reads(samp, metrics):
    numreads = metrics[metrics['cell_id'] == samp]['total_mapped_reads'].iloc[0]
    return numreads

def get_cell_quality(samp, metrics):
    quality = metrics[metrics['cell_id'] == samp]['probability_good'].iloc[0]
    return quality

def merge_pdf(in_filenames, out_filename, metrics, mad_threshold, numreads_threshold, quality_theshold):

    metrics = pd.read_csv(metrics, sep=',')
    expconds = metrics["experimental_condition"].unique()

    cells = []
    for expcond in expconds:
        cells += sorted(metrics[metrics["experimental_condition"]==expcond]["cell_id"])

    for infiles, out_file in zip(in_filenames, out_filename):

        helpers.makedirs(out_file, isfile=True)

        if mad_threshold:
            cells = [cell for cell in cells if get_mad_score(cell, metrics) <= mad_threshold]

        if numreads_threshold:
            cells = [cell for cell in cells if get_num_reads(cell, metrics) >= numreads_threshold]

        if quality_theshold:
            cells = [cell for cell in cells if get_cell_quality(cell, metrics) >= quality_theshold]

        infiles = [infiles[samp] for samp in cells]

        pdfutils.merge_pdfs(infiles, out_file)
