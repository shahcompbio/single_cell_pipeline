'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner
import pandas as pd
from scripts import ExtractHmmMetrics
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


scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
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


def run_hmmcopy_script(corrected_reads, hmm_reads_out, hmm_segs_out,
                       hmm_params_out, hmm_posteriors_out, cell_id, hmmparams):

    # run hmmcopy
    cmd = ['Rscript', run_hmmcopy_rscript,
           '--corrected_data=' + corrected_reads,
           '--reads_output=' + hmm_reads_out,
           '--segs_output=' + hmm_segs_out,
           '--params_output=' + hmm_params_out,
           '--post_marginals_output=' + hmm_posteriors_out,
           '--sample_id=' + cell_id]

    if hmmparams['map_cutoff']:
        cmd.append('--map_cutoff=' + str(hmmparams['map_cutoff']))

    if hmmparams['num_states']:
        cmd.append('--num_states=' + str(hmmparams['num_states']))

    if hmmparams['mu']:
        cmd.append('--param_mu=' + str(hmmparams['mu']))

    if hmmparams['m']:
        cmd.append('--param_m=' + str(hmmparams['m']))

    if hmmparams['kappa']:
        cmd.append('--param_k=' + str(hmmparams['kappa']))

    if hmmparams['e']:
        cmd.append('--param_e=' + str(hmmparams['e']))

    if hmmparams['g']:
        cmd.append('--param_g=' + str(hmmparams['g']))

    if hmmparams['s']:
        cmd.append('--param_s=' + str(hmmparams['s']))

    if hmmparams['strength']:
        cmd.append('--param_str=' + str(hmmparams['strength']))

    if hmmparams['nu']:
        cmd.append('--param_nu=' + str(hmmparams['nu']))

    if hmmparams['eta']:
        cmd.append('--param_eta=' + str(hmmparams['eta']))

    if hmmparams['lambda']:
        cmd.append('--param_l=' + str(hmmparams['lambda']))

    if hmmparams["auto_ploidy"]:
        cmd.append('--auto_ploidy')

    pypeliner.commandline.execute(*cmd)


def run_hmmcopy(
        bam_file,
        bai_file,
        corrected_reads_filename,
        segments_filename,
        parameters_filename,
        posterior_marginals_filename,
        hmm_metrics,
        reads_plot,
        segs_plot,
        bias_plot,
        params_plot,
        cell_id,
        config,
        hmmparams,
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

    hmm_reads_out = os.path.join(tempdir, "hmm_reads.csv")
    hmm_params_out = os.path.join(tempdir, "hmm_params.csv")
    hmm_segs_out = os.path.join(tempdir, "hmm_segs.csv")
    hmm_posteriors_out = os.path.join(tempdir, "hmm_posteriors.csv")

    run_hmmcopy_script(
        corrected_reads,
        hmm_reads_out,
        hmm_segs_out,
        hmm_params_out,
        hmm_posteriors_out,
        cell_id,
        hmmparams)

    # generate the metrics file for hmmcopy
    metrics = ExtractHmmMetrics(hmm_params_out, hmm_reads_out,
                                hmm_segs_out, hmm_metrics, cell_id,
                                table_name='hmmcopy/metrics/{}'.format(cell_id))
    metrics.main()

    # convert hmmcopy outputs to h5
    hdfutils.convert_csv_to_hdf(
        hmm_params_out,
        parameters_filename,
        "/hmmcopy/params/{}".format(cell_id))
    hdfutils.convert_csv_to_hdf(
        hmm_posteriors_out,
        posterior_marginals_filename,
        "/hmmcopy/posteriors/{}".format(cell_id))
    hdfutils.convert_csv_to_hdf(
        hmm_segs_out,
        segments_filename,
        "/hmmcopy/segments/{}".format(cell_id))
    hdfutils.convert_csv_to_hdf(
        hmm_reads_out,
        corrected_reads_filename,
        "/hmmcopy/reads/{}".format(cell_id))

    # plot
    annotation_cols = ['cell_call', 'experimental_condition', 'sample_type',
                       'coverage_depth', 'mad_neutral_state',
                       'MSRSI_non_integerness']
    plot_hmmcopy(
        corrected_reads_filename, segments_filename, parameters_filename,
        hmm_metrics, None, config["ref_genome"], reads_plot, segs_plot,
        bias_plot, params_plot, cell_id, hmmparams['num_states'], "QC", None,
        annotation_cols)


def annotate_metrics(metrics, sample_info, output):
    """
    adds sample information to metrics in place
    """
    with pd.HDFStore(output, 'w', complevel=9, complib='blosc') as output, pd.HDFStore(metrics) as metrics:
        for tableid in metrics.keys():
            data = metrics[tableid]

            cell_id = data["cell_id"].iloc[0]

            cell_info = sample_info[cell_id]

            for colname, value in cell_info.iteritems():
                data[colname] = value

                output.put(tableid, data, format="table")


def merge_files(reads, segs, hmm_metrics, hmm_params, hmm_posteriors,
                merged_segs, merged_reads, merged_hmm_metrics,
                merged_hmm_params, merged_hmm_posteriors,
                mad_thres, config, igv_segs, sample_info, tempdir):

    helpers.makedirs(tempdir)
    temp_hmm_metrics = os.path.join(tempdir, "merged_hmm_metrics.csv")

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
    hdfutils.concat_hdf_tables(
        hmm_posteriors,
        merged_hmm_posteriors,
        non_numeric_as_category=False)

    annotate_metrics(temp_hmm_metrics, sample_info, merged_hmm_metrics)

    convert_csv_to_seg(merged_segs, config["bin_size"], merged_hmm_metrics, igv_segs,
                       0.2)


def convert_csv_to_seg(segs, bin_size, metrics, output_seg, mad_threshold):
    converter = ConvertCSVToSEG(
        segs,
        bin_size,
        metrics,
        output_seg,
        mad_threshold)
    converter.main()


def plot_hmmcopy(reads, segments, params, metrics, sample_info, ref_genome, reads_out, segs_out,
                 bias_out, params_out, cell_id, num_states=7, plot_title=None,
                 mad_threshold=None, annotation_cols=None):
    plot = GenHmmPlots(reads, segments, params, metrics, sample_info, ref_genome, reads_out, segs_out,
                       bias_out, params_out, cell_id, num_states=num_states, plot_title=plot_title,
                       mad_threshold=mad_threshold, annotation_cols=annotation_cols)
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


def get_hierarchical_clustering_order(reads_filename, cluster_order_output):

    with pd.HDFStore(reads_filename, 'r') as reads_store:

        data = []

        for tableid in reads_store.keys():

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
    order = {v: i for i, v in enumerate(order)}

    return order


def add_clustering_order(
        reads_filename, metrics_filename, cluster_order_output):
    order = get_hierarchical_clustering_order(
        reads_filename,
        cluster_order_output)

    with pd.HDFStore(metrics_filename, 'r') as metrics, pd.HDFStore(cluster_order_output, 'w', complevel=9, complib='blosc') as output:

        for tableid in metrics.keys():
            data = metrics[tableid]

            cell_id = data["cell_id"].iloc[0]

            data["order"] = order[cell_id]

            output.put(tableid, data, format="table")


def plot_kernel_density(infile, output, sep, colname, plot_title):
    plot = PlotKernelDensity(infile, output, sep, colname, plot_title)
    plot.main()


def plot_metrics(metrics, output, plot_title,):
    plot = PlotMetrics(metrics, output, plot_title)
    plot.main()


def plot_pcolor(infile, metrics, output, plot_title=None,
                column_name=None, plot_by_col=None, numreads_threshold=None,
                mad_threshold=None, chromosomes=None, max_cn=None, median_hmmcopy_reads_per_bin_threshold=None):

    plot = PlotPcolor(infile, metrics, output, plot_title=plot_title,
                      column_name=column_name, plot_by_col=plot_by_col,
                      numreads_threshold=numreads_threshold,
                      mad_threshold=mad_threshold, chromosomes=chromosomes,
                      max_cn=max_cn,
                      median_hmmcopy_reads_per_bin_threshold=median_hmmcopy_reads_per_bin_threshold)
    plot.main()


def merge_cells_in_memory(hdf_input):
    data = []

    with pd.HDFStore(hdf_input, 'r') as input_store:
        for tableid in input_store:
            data.append(input_store[tableid])

    data = pd.concat(data)
    data = data.reset_index()

    return data


def merge_cells_on_disk(hdf_input, output_store_obj, tablename, dtypes={}):

    with pd.HDFStore(hdf_input, 'r') as input_store:
        for tableid in input_store:
            celldata = input_store[tableid]

            for col, dtype in dtypes.iteritems():
                celldata[col] = celldata[col].astype(dtype)

            if tablename not in output_store_obj:
                output_store_obj.put(tablename, celldata, format='table')
            else:
                output_store_obj.append(tablename, celldata, format='table')


def merge_stores(reads, segments, metrics, params, posteriors, output):

    with pd.HDFStore(output, 'w', complevel=9, complib='blosc') as output_store:
        output_store.put(
            '/hmmcopy/metrics',
            merge_cells_in_memory(metrics),
            format='table')
        output_store.put(
            '/hmmcopy/params',
            merge_cells_in_memory(params),
            format='table')
        output_store.put(
            '/hmmcopy/posteriors',
            merge_cells_in_memory(posteriors),
            format='table')

        dtypes = {'integer_copy_number': float, 'valid': bool, 'ideal': bool}
        merge_cells_on_disk(
            reads,
            output_store,
            '/hmmcopy/reads',
            dtypes=dtypes)
        dtypes = {'integer_copy_number': float}
        merge_cells_on_disk(
            segments,
            output_store,
            '/hmmcopy/segments',
            dtypes=dtypes)
