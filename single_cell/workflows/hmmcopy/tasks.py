'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import shutil
import pypeliner
import pandas as pd
import numpy as np
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc

from scripts import ConvertCSVToSEG
from scripts import ReadCounter
from scripts import CorrectReadCount
from scripts import classify

from single_cell.utils import pdfutils
from single_cell.utils import helpers
from single_cell.utils import hdfutils
from single_cell.utils import csvutils

from single_cell.utils.singlecell_copynumber_plot_utils import PlotKernelDensity
from single_cell.utils.singlecell_copynumber_plot_utils import PlotMetrics
from single_cell.utils.singlecell_copynumber_plot_utils import PlotPcolor
from single_cell.utils.singlecell_copynumber_plot_utils import GenHmmPlots


scripts_directory = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)),
    'scripts')
run_hmmcopy_rscript = os.path.join(scripts_directory, 'hmmcopy.R')


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
    cmd = ["hmmcopy"]

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
        segs_pdf_filename,
        bias_pdf_filename,
        cell_id,
        reference,
        hmmparams,
        multipliers,
        tempdir,
        sample_info,
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

    run_hmmcopy_script(
        corrected_reads,
        tempdir,
        cell_id,
        hmmparams,
        docker_image
    )

    for multiplier in multipliers:
        hmmcopy_outdir = os.path.join(tempdir, str(multiplier))

        hmmcopy_reads_file = os.path.join(hmmcopy_outdir, "reads.csv")
        hmmcopy_params_files = os.path.join(hmmcopy_outdir, "params.csv")
        hmmcopy_segs_files = os.path.join(hmmcopy_outdir, "segs.csv")
        hmmcopy_metrics_files = os.path.join(hmmcopy_outdir, "metrics.csv")

        gzip_file(hmmcopy_reads_file, corrected_reads_filename[multiplier])
        gzip_file(hmmcopy_params_files, parameters_filename[multiplier])
        gzip_file(hmmcopy_segs_files, segments_filename[multiplier])
        gzip_file(hmmcopy_metrics_files, metrics_filename[multiplier])


    corrected_reads_all = [corrected_reads_filename[mult] for mult in multipliers]
    segments_all = [segments_filename[mult] for mult in multipliers]
    params_all = [parameters_filename[mult] for mult in multipliers]
    metrics_all = [metrics_filename[mult] for mult in multipliers]

    annotation_cols = ['cell_call', 'experimental_condition', 'sample_type',
                       'mad_neutral_state', 'MSRSI_non_integerness',
                       'total_mapped_reads_hmmcopy']

    plot_hmmcopy(
        corrected_reads_all, segments_all, params_all,
        metrics_all, reference, segs_pdf_filename,
        bias_pdf_filename, cell_id, multipliers,
        num_states=hmmparams['num_states'],
        annotation_cols=annotation_cols,
        sample_info=sample_info)



def concatenate_csv(inputs, output, yamloutput, multipliers, cells, low_memory=False, quick=False):
    for multiplier in multipliers:
        mult_inputs = [inputs[(cell, multiplier)]for cell in cells]
        if yamloutput:
            yamlfile = yamloutput[multiplier]
        else:
            yamlfile = None
        if low_memory:
            if quick:
                csvutils.concatenate_csv_lowmem_quick(
                    mult_inputs, output[multiplier], yamlfile=yamlfile
                )
            else:
                csvutils.concatenate_csv_lowmem(
                    mult_inputs, output[multiplier], yamlfile=yamlfile
                )
        else:
            csvutils.concatenate_csv(
                mult_inputs, output[multiplier], yamlfile=yamlfile
            )


def concatenate_csv_single(inputs, output, yamloutput, multiplier, cells, low_memory=False, quick=False):
        mult_inputs = [inputs[(cell, multiplier)]for cell in cells]
        if low_memory:
            if quick:
                csvutils.concatenate_csv_lowmem_quick(
                    mult_inputs, output, yamlfile=yamloutput
                )
            else:
                csvutils.concatenate_csv_lowmem(
                    mult_inputs, output, yamlfile=yamloutput
                )
        else:
            csvutils.concatenate_csv(
                mult_inputs, output, yamlfile=yamloutput
            )



def annotate_metrics(
        reads, metrics, output, sample_info, cells, chromosomes=None, yamlfile=None):
    """
    adds sample information to metrics in place
    """

    metrics = pd.read_csv(metrics, compression='gzip')

    order = get_hierarchical_clustering_order(
        reads, chromosomes=chromosomes
    )

    for cellid, value in order.iteritems():
        cellinfo = sample_info[cellid]
        value.update(cellinfo)

        for cellid in cells:
            cell_info = order[cellid]
            for colname, value in cell_info.iteritems():
                metrics.loc[metrics["cell_id"] == cellid, colname] = value

    metrics.to_csv(output, na_rep='NA', index=False, compression='gzip')

    if yamlfile:
        csvutils.generate_dtype_yaml(metrics, yamlfile)


def merge_hdf_files_on_disk(
        reads, merged_reads, multipliers, tableprefix, dtypes={}):

    output_store = pd.HDFStore(merged_reads, 'w', complevel=9, complib='blosc')

    cells = reads.keys()

    min_itemsize = {'cell_id': max(map(len, cells))+2}

    for cellid, infile in reads.iteritems():
        with pd.HDFStore(infile, 'r') as infilestore:
            for multiplier in multipliers:
                tablename = '/{}/{}/{}'.format(tableprefix, cellid, multiplier)
                data = infilestore[tablename]

                data = data.dropna(axis=0, how='all')
                if data.empty:
                    continue

                # remove later, only added to fix a run
                if "cell_id" not in data:
                    data.cell_id = cellid

                # make cellid categorical
                data.cell_id = pd.Categorical(data.cell_id, cells)

                for col, dtype in dtypes.iteritems():
                    data[col] = data[col].astype(dtype)

                out_tablename = '/{}/{}'.format(tableprefix, multiplier)
                if out_tablename not in output_store:
                    output_store.put(
                        out_tablename,
                        data,
                        format='table',)
                else:
                    output_store.append(
                        out_tablename,
                        data,
                        format='table',)

    output_store.close()


def merge_hdf_files_in_memory(
        reads, merged_reads, multipliers, tableprefix, dtypes={}):

    output_store = pd.HDFStore(merged_reads, 'w', complevel=9, complib='blosc')

    cells = reads.keys()

    for multiplier in multipliers:
        all_cells_data = []
        for cellid, infile in reads.iteritems():

            tablename = '/{}/{}/{}'.format(tableprefix, cellid, multiplier)

            with pd.HDFStore(infile, 'r') as infilestore:
                data = infilestore[tablename]

            for col, dtype in dtypes.iteritems():
                data[col] = data[col].astype(dtype)

            all_cells_data.append(data)

        all_cells_data = pd.concat(all_cells_data)

        all_cells_data = all_cells_data.reset_index()

        out_tablename = '/{}/{}'.format(tableprefix, multiplier)

        # make cellid categorical
        all_cells_data.cell_id = pd.Categorical(all_cells_data.cell_id, cells)

        output_store.put(
            out_tablename,
            all_cells_data,
            format='table')

    output_store.close()


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
        for row, coldata in grouped_data.iteritems():
            coldata = sorted(coldata, key=lambda tup: tup[1])
            coldata = [cell for cell, col in coldata]
            grouped_data[row] = coldata

    return grouped_data


def merge_pdf(in_filenames, outfilenames, metrics, cell_filters, tempdir, labels):

    helpers.makedirs(tempdir)

    metrics = metrics[0]

    good_cells = get_good_cells(
        metrics, cell_filters
    )

    grouped_data = group_cells_by_row(
        good_cells, metrics, sort_by_col=True
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
                   igv_segs, config, multiplier):

    converter = ConvertCSVToSEG(
        merged_segs,
        config['bin_size'],
        merged_hmm_metrics,
        igv_segs,
        config['igv_segs_quality_threshold'], multiplier)
    converter.main()


def plot_hmmcopy(reads, segments, params, metrics, ref_genome, segs_out,
                 bias_out, cell_id, multiplier, num_states=7,
                 annotation_cols=None, sample_info=None):

    with GenHmmPlots(reads, segments, params, metrics, ref_genome, segs_out,
                     bias_out, cell_id, multiplier, num_states=num_states,
                     annotation_cols=annotation_cols,
                     sample_info=sample_info) as plot:
        plot.main()


def extract_cell_by_col(df, colname, colvalue, rowname):
    return df[df[colname] == colvalue][rowname].iloc[0]


def get_good_cells(metrics, cell_filters):

    metrics_data = pd.read_csv(metrics, compression='gzip')

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


def sort_cells(metrics, good_cells, tableid):

    with pd.HDFStore(metrics, 'r') as metrics:

        data = metrics[tableid]

    cells_order = {
        order: cell for cell,
        order in zip(
            data.cell_id,
            data.order)}

    ordervals = sorted(cells_order.keys())

    sorted_cells = [cells_order[ordval] for ordval in ordervals]

    sorted_cells = [cell for cell in sorted_cells if cell in good_cells]

    return sorted_cells


def sort_bins(bins, chromosomes):

    bins = bins.drop_duplicates()

    if not chromosomes:
        chromosomes = map(str, range(1, 23)) + ['X', 'Y']

    bins["chr"] = pd.Categorical(bins["chr"], chromosomes)

    bins = bins.sort_values(['start', ])

    bins = [tuple(v) for v in bins.values.tolist()]

    return bins


def get_hierarchical_clustering_order(
        reads_filename, chromosomes=None):

    data = []
    chunksize = 10 ** 5
    for chunk in pd.read_csv(
            reads_filename, chunksize=chunksize, compression='gzip'):

        chunk["bin"] = list(zip(chunk.chr, chunk.start, chunk.end))

        print len(set(zip(chunk.cell_id.tolist(), chunk.bin.tolist())))

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


def plot_metrics(metrics, output, plot_title, multiplier):

    mult_plot_title = '{}({})'.format(plot_title, multiplier)

    tablename = '/hmmcopy/metrics/{}'.format(multiplier)

    plot = PlotMetrics(
        metrics,
        output,
        mult_plot_title,
        tablename=tablename)
    plot.plot_hmmcopy_metrics()


def plot_kernel_density(
        infile, output, sep, colname, plot_title, multiplier):

    mult_plot_title = '{}({})'.format(plot_title, multiplier)

    plot = PlotKernelDensity(
        infile,
        output,
        sep,
        colname,
        mult_plot_title)
    plot.main()


def plot_pcolor(infile, metrics, output, multiplier, plot_title=None,
                column_name=None, plot_by_col=None,
                chromosomes=None, max_cn=None,
                scale_by_cells=None, color_by_col=None,
                cell_filters=None, mappability_threshold=None):

    cells = get_good_cells(metrics, cell_filters)

    mult_plot_title = '{}({})'.format(plot_title, multiplier)

    reads_tablename = '/hmmcopy/reads/{}'.format(multiplier)
    metrics_tablename = '/hmmcopy/metrics/{}'.format(multiplier)

    plot = PlotPcolor(infile, metrics, output, plot_title=mult_plot_title,
                      column_name=column_name, plot_by_col=plot_by_col,
                      chromosomes=chromosomes,
                      max_cn=max_cn, multiplier=multiplier,
                      scale_by_cells=scale_by_cells,
                      segs_tablename=reads_tablename,
                      metrics_tablename=metrics_tablename,
                      color_by_col=color_by_col,
                      cells=cells,
                      mappability_threshold=mappability_threshold)
    plot.main()


def merge_tables(reads, segments, metrics, params, output, cells):
    categories = {'cell_id': cells}
    hdfutils.concat_hdf_tables([reads, segments, metrics, params], output, categories)


def add_quality(hmmcopy_metrics, alignment_metrics, multipliers, output, training_data):

    hmmcopy_tables = ['/hmmcopy/metrics/{}'.format(mult) for mult in multipliers]

    model = classify.train_classifier(training_data)

    feature_names = model.feature_names_

    data = classify.load_data(hmmcopy_metrics, alignment_metrics,
                              hmmcopy_tables, '/alignment/metrics',
                              feature_names)

    for hmmcopy_table, tabledata in data:
        predictions = classify.classify(model, tabledata)

        classify.write_to_output(
            hmmcopy_metrics,
            hmmcopy_table,
            output,
            predictions)


