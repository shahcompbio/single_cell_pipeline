'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import shutil

import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hc
import scipy.spatial as sp
from single_cell.utils import csvutils
from single_cell.utils import helpers
from single_cell.utils.singlecell_copynumber_plot_utils import GenHmmPlots
from single_cell.utils.singlecell_copynumber_plot_utils import PlotKernelDensity
from single_cell.utils.singlecell_copynumber_plot_utils import PlotMetrics
from single_cell.utils.singlecell_copynumber_plot_utils import PlotPcolor
from single_cell.workflows.hmmcopy.dtypes import dtypes

import pypeliner
from .scripts import ConvertCSVToSEG
from .scripts import CorrectReadCount
from .scripts import ReadCounter
from .scripts import classify

scripts_directory = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)),
    'scripts')
run_hmmcopy_rscript = os.path.join(scripts_directory, 'hmmcopy.R')


def get_max_cn(reads):
    df = csvutils.read_csv_and_yaml(reads)
    max_cn = np.nanpercentile(df['copy'], 99)
    return max_cn


def run_correction_hmmcopy(
        bam_file, correct_reads_out, readcount_wig, hmmparams, docker_image):
    run_readcount_rscript = os.path.join(
        scripts_directory,
        'correct_read_count.R')

    rc = ReadCounter(bam_file, readcount_wig, hmmparams['bin_size'], hmmparams['chromosomes'],
                     hmmparams['min_mqual'], excluded=hmmparams['exclude_list'])
    rc.main()

    if hmmparams["smoothing_function"] == 'loess':
        cmd = ['Rscript', run_readcount_rscript,
               readcount_wig,
               hmmparams['gc_wig_file'],
               hmmparams['map_wig_file'],
               correct_reads_out
               ]
        pypeliner.commandline.execute(*cmd, docker_image=docker_image)
    elif hmmparams["smoothing_function"] == 'modal':
        CorrectReadCount(hmmparams["gc_wig_file"],
                         hmmparams['map_wig_file'],
                         readcount_wig,
                         correct_reads_out,
                         mappability=hmmparams['map_cutoff']).main()
    else:
        raise Exception(
            "smoothing function %s not supported. pipeline supports loess and modal" %
            hmmparams["smoothing_function"])

    return correct_reads_out


def run_hmmcopy_script(corrected_reads, tempdir, cell_id, hmmparams, docker_image):
    cmd = ["hmmcopy_single_cell.R"]

    # run hmmcopy
    cmd += ['--corrected_data=' + corrected_reads,
            '--outdir=' + tempdir,
            '--sample_id=' + cell_id]

    multipliers = ','.join(map(str, hmmparams['multipliers']))

    cmd.append('--param_str=' + str(hmmparams['strength']))
    cmd.append('--param_e=' + str(hmmparams['e']))
    cmd.append('--param_mu=' + str(hmmparams['mu']))
    cmd.append('--param_l=' + str(hmmparams['lambda']))
    cmd.append('--param_nu=' + str(hmmparams['nu']))
    cmd.append('--param_k=' + str(hmmparams['kappa']))
    cmd.append('--param_m=' + str(hmmparams['m']))
    cmd.append('--param_eta=' + str(hmmparams['eta']))
    cmd.append('--param_g=' + str(hmmparams['g']))
    cmd.append('--param_s=' + str(hmmparams['s']))
    cmd.append('--param_multiplier=' + multipliers)
    
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def gzip_file(inputfile, gzipped_csv):
    import gzip
    with open(inputfile) as inputdata, gzip.open(gzipped_csv, 'w') as gzipped_out:
        for line in inputdata:
            gzipped_out.write(line)


def run_hmmcopy(
        bam_file,
        corrected_reads_filename,
        segments_filename,
        parameters_filename,
        metrics_filename,
        hmmcopy_tar,
        cell_id,
        hmmparams,
        tempdir,
        docker_image
):

    # generate wig file for hmmcopy
    helpers.makedirs(tempdir)
    readcount_wig = os.path.join(tempdir, 'readcounter.wig')
    corrected_reads = os.path.join(tempdir, 'corrected_reads.csv')

    run_correction_hmmcopy(
        bam_file,
        corrected_reads,
        readcount_wig,
        hmmparams,
        docker_image
    )

    hmmcopy_tempdir = os.path.join(tempdir, '{}_hmmcopy'.format(cell_id))
    helpers.makedirs(hmmcopy_tempdir)

    run_hmmcopy_script(
        corrected_reads,
        hmmcopy_tempdir,
        cell_id,
        hmmparams,
        docker_image
    )

    hmmcopy_outdir = os.path.join(hmmcopy_tempdir, str(0))
    
    csvutils.rewrite_csv_file(
        os.path.join(hmmcopy_outdir, "reads.csv"), corrected_reads_filename,
        dtypes=dtypes()['reads']
    )
    
    csvutils.rewrite_csv_file(
        os.path.join(hmmcopy_outdir, "params.csv"), parameters_filename,
        dtypes=dtypes()['params']
    )
 
    csvutils.rewrite_csv_file(
        os.path.join(hmmcopy_outdir, "segs.csv"), segments_filename,
        dtypes=dtypes()['segs']
    )
    
    csvutils.rewrite_csv_file(
        os.path.join(hmmcopy_outdir, "metrics.csv"), metrics_filename,
        dtypes=dtypes()['metrics']
    )

    helpers.make_tarfile(hmmcopy_tar, hmmcopy_tempdir)


def concatenate_csv(inputs, output):
    csvutils.concatenate_csv(
        inputs,
        output,
        write_header=True
    )


def get_hierarchical_clustering_order(
        reads_filename, chromosomes=None):
    data = []
    chunksize = 10 ** 5
    for chunk in csvutils.read_csv_and_yaml(
            reads_filename, chunksize=chunksize):
        chunk["bin"] = list(zip(chunk.chr, chunk.start, chunk.end))

        # for some reason pivot doesnt like an Int64 state col
        chunk['state'] = chunk['state'].astype('float')

        chunk = chunk.pivot(index='cell_id', columns='bin', values='state')

        data.append(chunk)

    # merge chunks, sum cells that get split across chunks
    table = pd.concat(data)
    table = table.groupby(table.index).sum()

    bins = pd.DataFrame(
        table.columns.values.tolist(),
        columns=[
            'chr',
            'start',
            'end'])

    bins['chr'] = bins['chr'].astype(str)

    bins = sort_bins(bins, chromosomes)

    table = table.sort_values(bins, axis=0)

    data_mat = np.array(table.values)

    data_mat[np.isnan(data_mat)] = -1

    row_linkage = hc.linkage(sp.distance.pdist(data_mat, 'cityblock'),
                             method='ward')

    order = hc.leaves_list(row_linkage)

    samps = table.index
    order = [samps[i] for i in order]
    order = {v: i for i, v in enumerate(order)}

    return order


def get_mappability_col(reads, annotated_reads):
    reads = csvutils.read_csv_and_yaml(reads, chunksize=100)

    alldata = []
    for read_data in reads:
        read_data['is_low_mappability'] = (read_data['map'] <= 0.9)
        alldata.append(read_data)

    alldata = pd.concat(alldata)

    csvutils.write_dataframe_to_csv_and_yaml(
        alldata, annotated_reads, dtypes()['reads'], write_header=True
    )


def add_clustering_order(
        reads, metrics, output, chromosomes=None, sample_info=None):
    """
    adds sample information to metrics in place
    """

    order = get_hierarchical_clustering_order(
        reads, chromosomes=chromosomes
    )

    if not sample_info:
        sample_info = {}

    for cell_id, order in order.items():
        if cell_id not in sample_info:
            sample_info[cell_id] = {}
        sample_info[cell_id]['order'] = order

    csvutils.annotate_csv(metrics, sample_info, output, dtypes()['metrics'])

def group_cells_by_row(cells, metrics, sort_by_col=False):
    metricsdata = pd.read_csv(metrics, compression='gzip')
    metricsdata = metricsdata.set_index('cell_id')

    grouped_data = {}

    for cell in cells:

        row = metricsdata.loc[cell]["row"]
        col = metricsdata.loc[cell]["column"]

        if row not in grouped_data:
            grouped_data[row] = []

        value = (cell, col) if sort_by_col else cell

        grouped_data[row].append(value)

    if sort_by_col:
        for row, coldata in grouped_data.items():
            coldata = sorted(coldata, key=lambda tup: tup[1])
            coldata = [cell for cell, col in coldata]
            grouped_data[row] = coldata

    return grouped_data


def merge_pdf(in_filenames, outfilenames, metrics, cell_filters, tempdir, labels):
    helpers.makedirs(tempdir)

    good_cells = get_good_cells(
        metrics, cell_filters
    )

    for infiles, outfiles, label in zip(in_filenames, outfilenames, labels):

        extension = os.path.splitext(infiles[good_cells[0]])[-1]

        plotdir = os.path.join(tempdir, label)

        helpers.makedirs(plotdir)

        for cell in good_cells:
            shutil.copyfile(
                infiles[cell],
                os.path.join(plotdir, cell + "_" + label + extension))

        helpers.make_tarfile(outfiles, plotdir)


def create_igv_seg(merged_segs, merged_hmm_metrics,
                   igv_segs, config):
    converter = ConvertCSVToSEG(
        merged_segs,
        config['bin_size'],
        merged_hmm_metrics,
        igv_segs,
        config['igv_segs_quality_threshold'])
    converter.main()


def plot_hmmcopy(reads, segments, params, metrics, ref_genome, segs_out,
                 bias_out, cell_id, num_states=7,
                 annotation_cols=None, sample_info=None, max_cn=None):
    
    if not annotation_cols:
        annotation_cols = ['cell_call', 'experimental_condition', 'sample_type',
                           'mad_neutral_state', 'MSRSI_non_integerness',
                           'total_mapped_reads_hmmcopy']

    with GenHmmPlots(reads, segments, params, metrics, ref_genome, segs_out,
                     bias_out, cell_id, num_states=num_states,
                     annotation_cols=annotation_cols,
                     sample_info=sample_info, max_cn=max_cn) as plot:
        plot.main()


def get_good_cells(metrics, cell_filters):
    metrics_data = csvutils.read_csv_and_yaml(metrics)

    if not cell_filters:
        return metrics_data.cell_id.tolist()

    metrics_data = helpers.filter_metrics(metrics_data, cell_filters)

    return metrics_data.cell_id.tolist()


def sort_bins(bins, chromosomes):
    bins = bins.drop_duplicates()

    if not chromosomes:
        chromosomes = list(map(str, range(1, 23))) + ['X', 'Y']

    bins["chr"] = pd.Categorical(bins["chr"], chromosomes)

    bins = bins.sort_values(['start', ])

    bins = [tuple(v) for v in bins.values.tolist()]

    return bins


def plot_metrics(metrics, output, plot_title):
    plot = PlotMetrics(
        metrics,
        output,
        plot_title,
    )
    plot.plot_hmmcopy_metrics()


def plot_kernel_density(
        infile, output, sep, colname, plot_title):
    plot = PlotKernelDensity(
        infile,
        output,
        sep,
        colname,
        plot_title)
    plot.main()


def plot_pcolor(infile, metrics, output, plot_title=None,
                column_name=None, plot_by_col=None,
                chromosomes=None, max_cn=None,
                scale_by_cells=None, color_by_col=None,
                cell_filters=None, mappability_threshold=None):
    cells = get_good_cells(metrics, cell_filters)

    plot = PlotPcolor(
        infile, metrics, output, plot_title=plot_title,
        column_name=column_name, plot_by_col=plot_by_col,
        chromosomes=chromosomes,
        max_cn=max_cn,
        scale_by_cells=scale_by_cells,
        color_by_col=color_by_col,
        cells=cells,
        mappability_threshold=mappability_threshold
    )
    plot.main()


def add_quality(hmmcopy_metrics, alignment_metrics, multipliers, output, training_data, tempdir):
    helpers.makedirs(tempdir)

    hmmcopy_tables = ['/hmmcopy/metrics/{}'.format(mult) for mult in multipliers]

    model = classify.train_classifier(training_data)

    feature_names = model.feature_names_

    data = classify.load_data(hmmcopy_metrics, alignment_metrics,
                              hmmcopy_tables, '/alignment/metrics',
                              feature_names)

    for i, (hmmcopy_table, tabledata) in enumerate(data):
        intermediate_output = os.path.join(
            tempdir, '{}_metrics_with_quality.csv.gz'.format(i)
        )

        predictions = classify.classify(model, tabledata)

        classify.write_to_output(
            hmmcopy_metrics,
            hmmcopy_table,
            intermediate_output,
            predictions)

        csvutils.prep_csv_files(intermediate_output, output, dtypes=dtypes()['metrics'])
