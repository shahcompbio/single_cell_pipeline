'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner
from scripts import ReadCounter
from scripts import CorrectReadCount
from scripts import RunCopyClone

import pandas as pd
import numpy as np
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc

from single_cell.utils import csvutils
from single_cell.utils import pdfutils
from single_cell.utils import helpers
from single_cell.utils import hdfutils

from single_cell.utils.singlecell_copynumber_plot_utils import GenHmmPlots
from single_cell.utils.singlecell_copynumber_plot_utils import PlotMetrics
from single_cell.utils.singlecell_copynumber_plot_utils import PlotKernelDensity
from single_cell.utils.singlecell_copynumber_plot_utils import PlotPcolor


scripts_directory = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)),
    'scripts')
run_hmmcopy_rscript = os.path.join(scripts_directory, 'hmmcopy.R')


def merge_reads(reads, merged_reads):
    csvutils.concatenate_csv_lowmem(reads, merged_reads)


def correct_reads(
        bam_file,
        bai_file,
        reads_filename,
        config,
        tempdir,
        sample_id):

    ccparams = config["copyclone"]
    run_readcount_rscript = os.path.join(
        scripts_directory,
        'correct_read_count.R')

    # generate wig file for hmmcopy
    helpers.makedirs(tempdir)
    readcount_wig = os.path.join(tempdir, 'readcounter.wig')

    rc = ReadCounter(bam_file, readcount_wig, ccparams['bin_size'], config['chromosomes'],
                     ccparams['min_mqual'], excluded=ccparams['exclude_list'])
    rc.main()

    if ccparams["smoothing_function"] == 'loess':
        cmd = ['Rscript', run_readcount_rscript,
               readcount_wig,
               ccparams['gc_wig_file'],
               ccparams['map_wig_file'],
               reads_filename,
               ccparams["map_cutoff"],
               sample_id
               ]

        pypeliner.commandline.execute(*cmd)
    elif ccparams["smoothing_function"] == 'modal':
        CorrectReadCount(ccparams["gc_wig_file"],
                         ccparams['map_wig_file'],
                         readcount_wig,
                         reads_filename,
                         mappability=ccparams['map_cutoff'],
                         sample_id=sample_id).main()
    else:
        raise Exception(
            "smoothing function %s not supported. pipeline supports loess and modal" %
            ccparams["smoothing_function"])


def run_copyclone(corrected_data, reads, segments, metrics, segments_plot, bias_plot,
                  tempdir, config,
                  cells, sample_info, A=None, alpha_A=None,
                  tau=None, nu=None, eta=None, shape=None,
                  rate=None, ploidy_states=None, num_states=None, num_cores=1):

    helpers.makedirs(tempdir)

    csv_reads = os.path.join(tempdir, "reads.csv")
    csv_segments = os.path.join(tempdir, "segments.csv")
    csv_metrics = os.path.join(tempdir, "metrics.csv")

    rc = RunCopyClone(
        corrected_data, csv_reads, csv_segments, csv_metrics,
        tau=tau, nu=nu, eta=eta, shape=shape, rate=rate,
        ploidy_states=ploidy_states, num_states=num_states)
    rc.main()

    dtypes = {
        "ideal": float,
        "valid": float,
        "cell_id": "category",
        "chr": "category"}
    hdfutils.convert_csv_to_hdf(
        csv_reads,
        reads,
        '/copyclone/reads/',
        dtypes=dtypes)
    hdfutils.convert_csv_to_hdf(csv_segments, segments, '/copyclone/segments/')
    hdfutils.convert_csv_to_hdf(csv_metrics, metrics, '/copyclone/metrics/')

    annotation_cols = ['cell_call', 'experimental_condition', 'sample_type',
                       'mad_neutral_state', 'MSRSI_non_integerness',
                       'total_mapped_reads']

    args = []

    segsfiles = []
    biasfiles = []

    for cell in cells:
        segsout = os.path.join(tempdir, "{}_segs.pdf".format(cell))
        biasout = os.path.join(tempdir, "{}_bias.pdf".format(cell))
        arg = (reads, segments, metrics, config['ref_genome'], segsout, biasout,
               cell, num_states, annotation_cols, sample_info[cell])
        args.append(arg)

        segsfiles.append(segsout)
        biasfiles.append(biasout)

    helpers.run_in_parallel(_plot_hmmcopy_worker, args, ncores=num_cores)

    pdfutils.merge_pdfs(segsfiles, segments_plot)
    pdfutils.merge_pdfs(biasfiles, bias_plot)


def _plot_hmmcopy_worker(reads, segments, metrics, ref_genome, segsout,
                         biasout, cell, numstates, annotation_cols, sample_info):

    tablename_format = ['copyclone', 'type']

    with GenHmmPlots(reads, segments, None, metrics, ref_genome, segsout,
                     biasout, cell, None, tablename_format,
                     num_states=numstates,
                     annotation_cols=annotation_cols,
                     sample_info=sample_info) as plot:
        plot.main()


def plot_metrics(metrics, output, tempdir, plot_title):

    tablename = '/copyclone/metrics/'

    plot = PlotMetrics(
        metrics,
        output,
        plot_title,
        tablename=tablename)
    plot.plot_hmmcopy_metrics()


def annotate_metrics(
        reads, metrics, output, sample_info, cells, chromosomes=None):
    """
    adds sample information to metrics in place
    """

    metrics_store = pd.HDFStore(metrics, 'r')

    output_store = pd.HDFStore(output, 'w', complevel=9, complib='blosc')

    tablename = '/copyclone/reads/'
    order = get_hierarchical_clustering_order(
        reads, tablename, chromosomes=chromosomes
    )
    for cellid, value in order.iteritems():
        cellinfo = sample_info[cellid]
        value.update(cellinfo)

    tablename = '/copyclone/metrics/'

    data = metrics_store[tablename]

    for cellid in cells:
        cell_info = order[cellid]
        for colname, value in cell_info.iteritems():
            data.loc[data["cell_id"] == cellid, colname] = value

    output_store.put(tablename, data, format='table')

    output_store.close()
    metrics_store.close()


def get_hierarchical_clustering_order(
        reads_filename, tablename, chromosomes=None):

    print reads_filename

    data = []
    chunksize = 10 ** 6
    for chunk in pd.read_hdf(
            reads_filename, chunksize=chunksize, key=tablename):

        chunk["bin"] = list(zip(chunk.chr, chunk.start, chunk.end))

        chunk = chunk.pivot(index='cell_id', columns='bin', values='state')

        data.append(chunk)

    table = pd.concat(data)

    bins = pd.DataFrame(
        table.columns.values.tolist(),
        columns=[
            'chr',
            'start',
            'end'])

    bins = sort_bins(bins, chromosomes)

    table = table.sort_values(bins, axis=0)

    data_mat = np.array(table.values)

    data_mat[np.isnan(data_mat)] = -1

    row_linkage = hc.linkage(sp.distance.pdist(data_mat, 'cityblock'),
                             method='ward')

    order = hc.leaves_list(row_linkage)

    samps = table.index
    order = [samps[i] for i in order]
    order = {v: {"order": i} for i, v in enumerate(order)}

    return order


def sort_bins(bins, chromosomes):

    bins = bins.drop_duplicates()

    if not chromosomes:
        chromosomes = map(str, range(1, 23)) + ['X', 'Y']

    bins["chr"] = pd.Categorical(bins["chr"], chromosomes)

    bins = bins.sort_values(['start', ])

    bins = [tuple(v) for v in bins.values.tolist()]

    return bins


def plot_kernel_density(
        infile, output, sep, colname, plot_title):

    tablename = '/copyclone/metrics/'

    plot = PlotKernelDensity(
        infile,
        output,
        sep,
        colname,
        plot_title,
        tablename=tablename)
    plot.main()


def extract_cell_by_col(df, colname, colvalue, rowname):
    return df[df[colname] == colvalue][rowname].iloc[0]


def get_good_cells(metrics, cell_filters, tableid):

    with pd.HDFStore(metrics, 'r') as metrics_store:
        metrics_data = metrics_store[tableid]

    cells = metrics_data.cell_id

    if not cell_filters:
        return cells.tolist()

    # cells to keep
    for metric_col, operation, threshold in cell_filters:
        cells = [cell for cell in cells
                 if helpers.eval_expr(
                     extract_cell_by_col(
                         metrics_data,
                         'cell_id',
                         cell,
                         metric_col),
                     operation,
                     threshold)]

    return cells


def plot_pcolor(infile, metrics, output, plot_title=None,
                column_name=None, plot_by_col=None,
                chromosomes=None, max_cn=None,
                scale_by_cells=None, color_by_col=None,
                cell_filters=None, mappability_threshold=None):

    cells = get_good_cells(metrics, cell_filters, '/copyclone/metrics/')

    reads_tablename = '/copyclone/reads/'
    metrics_tablename = '/copyclone/metrics/'

    plot = PlotPcolor(infile, metrics, output, plot_title=plot_title,
                      column_name=column_name, plot_by_col=plot_by_col,
                      chromosomes=chromosomes,
                      max_cn=max_cn,
                      scale_by_cells=scale_by_cells,
                      segs_tablename=reads_tablename,
                      metrics_tablename=metrics_tablename,
                      color_by_col=color_by_col,
                      cells=cells,
                      mappability_threshold=mappability_threshold)
    plot.main()
