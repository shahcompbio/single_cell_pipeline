'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner
import pandas as pd
from scripts import ConvertCSVToSEG
from scripts import ReadCounter
from scripts import CorrectReadCount

from single_cell.utils import pdfutils
from single_cell.utils import helpers
from single_cell.utils import hdfutils


from single_cell.utils.singlecell_copynumber_plot_utils import PlotKernelDensity
from single_cell.utils.singlecell_copynumber_plot_utils import PlotMetrics
from single_cell.utils.singlecell_copynumber_plot_utils import PlotPcolor
from single_cell.utils.singlecell_copynumber_plot_utils import GenHmmPlots


import scipy.spatial as sp
import scipy.cluster.hierarchy as hc


scripts_directory = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)),
    'scripts')
run_hmmcopy_rscript = os.path.join(scripts_directory, 'hmmcopy.R')


def run_correction_hmmcopy(
        bam_file, correct_reads_out, readcount_wig, config, hmmparams):

    run_readcount_rscript = os.path.join(
        scripts_directory,
        'correct_read_count.R')

    rc = ReadCounter(bam_file, readcount_wig, hmmparams['bin_size'], config['chromosomes'],
                     hmmparams['min_mqual'], excluded=hmmparams['exclude_list'])
    rc.main()

    if hmmparams["smoothing_function"] == 'loess':
        cmd = ['Rscript', run_readcount_rscript,
               readcount_wig,
               hmmparams['gc_wig_file'],
               hmmparams['map_wig_file'],
               correct_reads_out
               ]
        pypeliner.commandline.execute(*cmd)
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


def run_hmmcopy_script(corrected_reads, tempdir, cell_id, hmmparams):

    # run hmmcopy
    cmd = ['Rscript', run_hmmcopy_rscript,
           '--corrected_data=' + corrected_reads,
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

    pypeliner.commandline.execute(*cmd)


def run_hmmcopy(
        bam_file,
        bai_file,
        corrected_reads_filename,
        segments_filename,
        parameters_filename,
        metrics_filename,
        segs_pdf_filename,
        bias_pdf_filename,
        cell_id,
        config,
        hmmparams,
        multipliers,
        tempdir,):

    # generate wig file for hmmcopy
    os.makedirs(tempdir)
    readcount_wig = os.path.join(tempdir, 'readcounter.wig')
    corrected_reads = os.path.join(tempdir, 'corrected_reads.csv')

    run_correction_hmmcopy(
        bam_file,
        corrected_reads,
        readcount_wig,
        config,
        hmmparams)

    run_hmmcopy_script(
        corrected_reads,
        tempdir,
        cell_id,
        hmmparams)

    hmmcopy_reads_files = []
    hmmcopy_params_files = []
    hmmcopy_segs_files = []
    hmmcopy_metrics_files = []

    hmmcopy_reads_tablenames = []
    hmmcopy_params_tablenames = []
    hmmcopy_segs_tablenames = []
    hmmcopy_metrics_tablenames = []

    for multiplier in multipliers:
        hmmcopy_outdir = os.path.join(tempdir, str(multiplier))

        hmmcopy_reads_files.append(os.path.join(hmmcopy_outdir, "reads.csv"))
        hmmcopy_params_files.append(os.path.join(hmmcopy_outdir, "params.csv"))
        hmmcopy_segs_files.append(os.path.join(hmmcopy_outdir, "segs.csv"))
        hmmcopy_metrics_files.append(
            os.path.join(
                hmmcopy_outdir,
                "metrics.csv"))

        hmmcopy_reads_tablenames.append(
            '/hmmcopy/reads/{}/{}'.format(cell_id, multiplier))
        hmmcopy_params_tablenames.append(
            '/hmmcopy/params/{}/{}'.format(cell_id, multiplier))
        hmmcopy_segs_tablenames.append(
            '/hmmcopy/segments/{}/{}'.format(cell_id, multiplier))
        hmmcopy_metrics_tablenames.append(
            '/hmmcopy/metrics/{}/{}'.format(cell_id, multiplier))

    hdfutils.concat_csvs_to_hdf(
        hmmcopy_reads_files,
        corrected_reads_filename,
        hmmcopy_reads_tablenames)
    hdfutils.concat_csvs_to_hdf(
        hmmcopy_params_files,
        parameters_filename,
        hmmcopy_params_tablenames)
    hdfutils.concat_csvs_to_hdf(
        hmmcopy_segs_files,
        segments_filename,
        hmmcopy_segs_tablenames)
    hdfutils.concat_csvs_to_hdf(
        hmmcopy_metrics_files,
        metrics_filename,
        hmmcopy_metrics_tablenames)

    annotation_cols = ['cell_call', 'experimental_condition', 'sample_type',
                       'coverage_depth', 'mad_neutral_state',
                       'MSRSI_non_integerness']

    plot_hmmcopy(
        corrected_reads_filename, segments_filename, parameters_filename,
        metrics_filename, config["ref_genome"], segs_pdf_filename,
        bias_pdf_filename, cell_id, multipliers,
        num_states=hmmparams['num_states'],
        annotation_cols=annotation_cols)


def annotate_metrics(metrics, sample_info, output):
    """
    adds sample information to metrics in place
    """
    hdfutils.annotate_per_cell_store_with_dict(metrics, sample_info, output)


def merge_files(reads, segs, hmm_metrics, hmm_params,
                merged_segs, merged_reads, merged_hmm_metrics,
                merged_hmm_params, tempdir, igv_segs, sample_info,
                mad_thres, config, multipliers):

    helpers.makedirs(tempdir)
    temp_hmm_metrics = os.path.join(tempdir, "merged_hmm_metrics.h5")

    hdfutils.concat_hdf_tables(
        reads,
        merged_reads,
        in_memory=False,
        non_numeric_as_category=False)
    hdfutils.concat_hdf_tables(
        segs,
        merged_segs,
        in_memory=False,
        non_numeric_as_category=False)
    hdfutils.concat_hdf_tables(
        hmm_metrics,
        temp_hmm_metrics,
        non_numeric_as_category=False)

    hdfutils.concat_hdf_tables(
        hmm_params,
        merged_hmm_params,
        non_numeric_as_category=False)

    annotate_metrics(temp_hmm_metrics, sample_info, merged_hmm_metrics)

    convert_csv_to_seg(merged_segs, config["bin_size"], merged_hmm_metrics, igv_segs,
                       0.2, 0)


def convert_csv_to_seg(
        segs, bin_size, metrics, output_seg, mad_threshold, multiplier):
    converter = ConvertCSVToSEG(
        segs,
        bin_size,
        metrics,
        output_seg,
        mad_threshold, multiplier)
    converter.main()


def plot_hmmcopy(reads, segments, params, metrics, ref_genome, segs_out,
                 bias_out, cell_id, multiplier, num_states=7,
                 annotation_cols=None):

    with GenHmmPlots(reads, segments, params, metrics, ref_genome, segs_out,
                     bias_out, cell_id, multiplier, num_states=num_states,
                     annotation_cols=annotation_cols) as plot:
        plot.main()


def get_good_cells(metrics, mad_threshold, numreads_threshold,
                   median_hmmcopy_reads_per_bin_threshold):

    good_cells = []

    with pd.HDFStore(metrics, 'r') as metrics:
        for tableid in metrics.keys():
            data = metrics[tableid]

            cellid = data["cell_id"].iloc[0]

            mad_score = data["mad_neutral_state"].iloc[0]

            if mad_score > mad_threshold:
                continue

            num_reads = data['total_mapped_reads'].iloc[0]

            if num_reads < numreads_threshold:
                continue

            reads_per_bin = data['median_hmmcopy_reads_per_bin'].iloc[0]
            if reads_per_bin < median_hmmcopy_reads_per_bin_threshold:
                continue

            good_cells.append(cellid)

    return good_cells


def sort_cells(metrics, good_cells):

    cells_by_condition = {}

    with pd.HDFStore(metrics, 'r') as metrics:
        for tableid in metrics.keys():
            data = metrics[tableid]

            cellid = data["cell_id"].iloc[0]
            expcond = data["experimental_condition"].iloc[0]

            if expcond not in cells_by_condition:
                cells_by_condition[expcond] = []

            cells_by_condition[expcond].append(cellid)

    sorted_cells = []

    for expcond, cells in cells_by_condition.iteritems():
        pass_cells = [cell for cell in cells if cell in good_cells]

        sorted_cells += sorted(pass_cells)

    return sorted_cells


def merge_pdf(in_filenames, out_filename, metrics, mad_threshold,
              numreads_threshold, median_hmmcopy_reads_per_bin_threshold):

    good_cells = get_good_cells(
        metrics,
        mad_threshold,
        numreads_threshold,
        median_hmmcopy_reads_per_bin_threshold)

    sorted_cells = sort_cells(metrics, good_cells)

    for infiles, out_file in zip(in_filenames, out_filename):

        helpers.makedirs(out_file, isfile=True)

        infiles = [infiles[samp] for samp in sorted_cells]

        pdfutils.merge_pdfs(infiles, out_file)


def sort_bins(bins):

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


def get_hierarchical_clustering_order(
        reads_filename, cluster_order_output, multiplier):

    with pd.HDFStore(reads_filename, 'r') as reads_store:

        data = []

        for tableid in reads_store.keys():

            table_multiplier = int(tableid.split('/')[-1])
            if not table_multiplier == multiplier:
                continue

            table = reads_store[tableid]

            bins = table[['chr', 'start', 'end']]

            bins = sort_bins(bins)

            table["bin"] = list(zip(table.chr, table.start, table.end))

            table = table.pivot(index='cell_id', columns='bin', values='state')

            table = table.sort_values(bins, axis=0)

            data.append(table)

    data = pd.concat(data)

    row_linkage = hc.linkage(sp.distance.pdist(data.values),
                             method='average')
    order = hc.leaves_list(row_linkage)

    samps = data.index
    order = [samps[i] for i in order]
    order = {v: {"order": i} for i, v in enumerate(order)}

    return order


def add_clustering_order(
        reads_filename, metrics_filename, cluster_order_output, multipliers, cells):

    with pd.HDFStore(cluster_order_output, 'w', complevel=9, complib='blosc') as out_store:

        for multiplier in multipliers:
            order = get_hierarchical_clustering_order(
                reads_filename,
                cluster_order_output, multiplier)

            tables = [
                '/hmmcopy/metrics/{}/{}'.format(cell, multiplier) for cell in cells]

            hdfutils.annotate_store_with_dict(
                metrics_filename,
                order,
                out_store,
                tables=tables)


def plot_kernel_density(
        infile, output, tempdir, sep, colname, plot_title, multipliers):

    helpers.makedirs(tempdir)

    multiplier_pdfs = []

    for multiplier in multipliers:

        multiplier_output = os.path.join(tempdir, "{}.pdf".format(multiplier))

        multiplier_pdfs.append(multiplier_output)

        mult_plot_title = '{}({})'.format(plot_title, multiplier)

        plot = PlotKernelDensity(
            infile,
            multiplier_output,
            sep,
            colname,
            mult_plot_title,
            multiplier=multiplier)
        plot.main()

    pdfutils.merge_pdfs(multiplier_pdfs, output)


def plot_pcolor(infile, metrics, output, tempdir, multipliers, plot_title=None,
                column_name=None, plot_by_col=None, numreads_threshold=None,
                mad_threshold=None, chromosomes=None, max_cn=None,
                median_hmmcopy_reads_per_bin_threshold=None,
                ):

    helpers.makedirs(tempdir)

    multiplier_pdfs = []

    for multiplier in multipliers:

        multiplier_output = os.path.join(tempdir, "{}.pdf".format(multiplier))

        multiplier_pdfs.append(multiplier_output)

        mult_plot_title = '{}({})'.format(plot_title, multiplier)

        plot = PlotPcolor(infile, metrics, multiplier_output, plot_title=mult_plot_title,
                          column_name=column_name, plot_by_col=plot_by_col,
                          numreads_threshold=numreads_threshold,
                          mad_threshold=mad_threshold, chromosomes=chromosomes,
                          max_cn=max_cn, multiplier=multiplier,
                          median_hmmcopy_reads_per_bin_threshold=median_hmmcopy_reads_per_bin_threshold)
        plot.main()

    pdfutils.merge_pdfs(multiplier_pdfs, output)


def merge_tables(reads, segments, metrics, params, output, multipliers, cells):

    with pd.HDFStore(output, 'w', complevel=9, complib='blosc') as output_store:

        for multiplier in multipliers:

            tables_to_merge = [
                '/hmmcopy/metrics/{}/{}'.format(cell, multiplier) for cell in cells]
            hdfutils.merge_per_cell_tables(
                metrics,
                output_store,
                '/hmmcopy/metrics/{}'.format(multiplier),
                in_memory=True,
                dtypes={'too_even': float},
                tables_to_merge=tables_to_merge)

            tables_to_merge = [
                '/hmmcopy/params/{}/{}'.format(cell, multiplier) for cell in cells]
            hdfutils.merge_per_cell_tables(
                params,
                output_store,
                '/hmmcopy/params/{}'.format(multiplier),
                in_memory=True,
                tables_to_merge=tables_to_merge)

            dtypes = {'valid': bool, 'ideal': bool}
            tables_to_merge = [
                '/hmmcopy/reads/{}/{}'.format(cell, multiplier) for cell in cells]
            hdfutils.merge_per_cell_tables(
                reads,
                output_store,
                '/hmmcopy/reads/{}'.format(multiplier),
                in_memory=False,
                dtypes=dtypes,
                tables_to_merge=tables_to_merge)
            dtypes = {}
            tables_to_merge = [
                '/hmmcopy/segments/{}/{}'.format(cell, multiplier) for cell in cells]
            hdfutils.merge_per_cell_tables(
                segments,
                output_store,
                '/hmmcopy/segments/{}'.format(multiplier),
                in_memory=False,
                dtypes=dtypes,
                tables_to_merge=tables_to_merge)


def plot_metrics(metrics, output, tempdir, plot_title, multipliers):

    helpers.makedirs(tempdir)

    multiplier_pdfs = []

    for multiplier in multipliers:

        multiplier_output = os.path.join(tempdir, "{}.pdf".format(multiplier))

        multiplier_pdfs.append(multiplier_output)

        mult_plot_title = '{}({})'.format(plot_title, multiplier)

        plot = PlotMetrics(
            metrics,
            multiplier_output,
            mult_plot_title,
            multiplier=multiplier)
        plot.plot_hmmcopy_metrics()

    pdfutils.merge_pdfs(multiplier_pdfs, output)
