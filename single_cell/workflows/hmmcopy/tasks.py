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


def run_hmmcopy_script(corrected_reads, hmm_reads_out, hmm_segs_out,
                       hmm_params_out, hmm_metrics_out, cell_id, hmmparams):

    # run hmmcopy
    cmd = ['Rscript', run_hmmcopy_rscript,
           '--corrected_data=' + corrected_reads,
           '--reads_output=' + hmm_reads_out,
           '--segs_output=' + hmm_segs_out,
           '--params_output=' + hmm_params_out,
           '--metrics_output=' + hmm_metrics_out,
           '--sample_id=' + cell_id]

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
    cmd.append('--param_multiplier='+str(hmmparams['multiplier']))

    pypeliner.commandline.execute(*cmd)


def run_hmmcopy(
        bam_file,
        bai_file,
        corrected_reads_filename,
        segments_filename,
        parameters_filename,
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
    hmm_metrics_out = os.path.join(tempdir, "hmm_metrics.csv")

    run_hmmcopy_script(
        corrected_reads,
        hmm_reads_out,
        hmm_segs_out,
        hmm_params_out,
        hmm_metrics_out,
        cell_id,
        hmmparams)

    # convert hmmcopy outputs to h5
    hdfutils.convert_csv_to_hdf(
        hmm_params_out,
        parameters_filename,
        "/hmmcopy/params/{}".format(cell_id))
    hdfutils.convert_csv_to_hdf(
        hmm_segs_out,
        segments_filename,
        "/hmmcopy/segments/{}".format(cell_id))
    hdfutils.convert_csv_to_hdf(
        hmm_reads_out,
        corrected_reads_filename,
        "/hmmcopy/reads/{}".format(cell_id))
    hdfutils.convert_csv_to_hdf(
        hmm_metrics_out,
        hmm_metrics,
        "/hmmcopy/metrics/{}".format(cell_id))

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
    hdfutils.annotate_per_cell_store_with_dict(metrics, sample_info, output)
    
def merge_files(reads, segs, hmm_metrics, hmm_params,
                merged_segs, merged_reads, merged_hmm_metrics,
                merged_hmm_params,
                mad_thres, config, igv_segs, sample_info, tempdir):

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
    order = {v: {"order":i} for i, v in enumerate(order)}

    return order


def add_clustering_order(
        reads_filename, metrics_filename, cluster_order_output):
    order = get_hierarchical_clustering_order(
        reads_filename,
        cluster_order_output)

    hdfutils.annotate_store_with_dict(metrics_filename, order, cluster_order_output)

def plot_kernel_density(infile, output, sep, colname, plot_title):
    plot = PlotKernelDensity(infile, output, sep, colname, plot_title)
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


def merge_tables(reads, segments, metrics, params, output):

    with pd.HDFStore(output, 'w', complevel=9, complib='blosc') as output_store:
        hdfutils.merge_per_cell_tables(
            metrics,
            output_store,
            '/hmmcopy/metrics',
            in_memory=True,
            dtypes={'too_even':float})
        hdfutils.merge_per_cell_tables(
            params,
            output_store,
            '/hmmcopy/params',
            in_memory=True)

        dtypes = {'valid': bool, 'ideal': bool}
        hdfutils.merge_per_cell_tables(
            reads,
            output_store,
            '/hmmcopy/reads',
            in_memory=False,
            dtypes=dtypes)
        dtypes = {}
        hdfutils.merge_per_cell_tables(
            segments,
            output_store,
            '/hmmcopy/segments',
            in_memory=False,
            dtypes=dtypes)


def plot_metrics(metrics, output, plot_title,):
    plot = PlotMetrics(metrics, output, plot_title)
    plot.plot_hmmcopy_metrics()

