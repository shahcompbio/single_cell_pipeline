import os
import shutil
import warnings
from scripts import CollectMetrics
from scripts import GenerateCNMatrix
from scripts import RunTrimGalore
from scripts import SummaryMetrics

from single_cell.utils import picardutils
from single_cell.utils import bamutils
from single_cell.utils import helpers
from single_cell.utils import csvutils
from single_cell.utils import gatkutils
from single_cell.utils import hdfutils

from single_cell.utils.singlecell_copynumber_plot_utils import PlotMetrics


def get_postprocess_metrics(infile, infile_bai, tempdir,
                    config, markdups_metrics, flagstat_metrics):

    if not os.path.exists(tempdir):
        helpers.makedirs(tempdir)

    outfile = os.path.join(tempdir, 'markdps.bam')
    outfile_index = outfile + '.bai'

    container_ctx = helpers.get_container_ctx(config['containers'], 'picard', docker_only=True)

    picardutils.bam_markdups(infile, outfile, markdups_metrics, tempdir,
                             **container_ctx)

    container_ctx = helpers.get_container_ctx(config['containers'], 'samtools', docker_only=True)
    bamutils.bam_index(outfile, outfile_index, **container_ctx)
    bamutils.bam_flagstat(outfile, flagstat_metrics, **container_ctx)

def plot_metrics(metrics, output, plot_title, gc_matrix, gc_content):

    plot = PlotMetrics(
        metrics, output, plot_title, gcbias_matrix=gc_matrix,
        gc_content=gc_content, tablename='/alignment/metrics',
        gc_tablename='/alignment/gc_metrics')
    plot.plot_alignment_metrics()


def get_summary_metrics(infile, output):
    summ = SummaryMetrics(infile, output)
    summ.main()


def annotate_metrics(infile, sample_info, outfile):

    hdfutils.annotate_store_with_dict(infile, sample_info, outfile)


def merge_all_metrics(infiles, outfile):
    csvutils.merge_csv(infiles, outfile, "outer", "cell_id")


def bam_collect_wgs_metrics(
        bam_filename, ref_genome, metrics_filename, config, tempdir):

    container_ctx = helpers.get_container_ctx(config['containers'], 'picard', docker_only=True)

    picardutils.bam_collect_wgs_metrics(
        bam_filename,
        ref_genome,
        metrics_filename,
        config,
        tempdir,
        **container_ctx)


def bam_collect_gc_metrics(
        bam_filename, ref_genome, metrics_filename, summary_filename, chart_filename, tempdir, config):
    container_ctx = helpers.get_container_ctx(config['containers'], 'picard', docker_only=True)

    picardutils.bam_collect_gc_metrics(
        bam_filename,
        ref_genome,
        metrics_filename,
        summary_filename,
        chart_filename,
        tempdir,
        **container_ctx)



def bam_collect_insert_metrics(
        bam_filename, flagstat_metrics_filename, metrics_filename, histogram_filename, tempdir,config):
    container_ctx = helpers.get_container_ctx(config['containers'], 'picard', docker_only=True)

    picardutils.bam_collect_insert_metrics(
        bam_filename,
        flagstat_metrics_filename,
        metrics_filename,
        histogram_filename,
        tempdir, **container_ctx)


def collect_gc(infiles, outfile, tempdir):

    helpers.makedirs(tempdir)

    tempouts = []
    for cell_id, infile in infiles.iteritems():
        tempout = os.path.join(
            tempdir,
            os.path.basename(infile) +
            ".parsed.csv")
        tempouts.append(tempout)
        gen_gc = GenerateCNMatrix(infile, tempout, ',',
                                  'NORMALIZED_COVERAGE', cell_id,
                                  'gcbias')
        gen_gc.main()

    merged_csv = os.path.join(tempdir, "merged_gc_metrics.csv")
    csvutils.concatenate_csv(tempouts, merged_csv)

    hdfutils.convert_csv_to_hdf(merged_csv, outfile, '/alignment/gc_metrics')


def collect_metrics(flagstat_metrics, markdups_metrics, insert_metrics,
                    wgs_metrics, tempdir, merged_metrics):

    helpers.makedirs(tempdir)
    sample_outputs = []
    for sample in flagstat_metrics.keys():
        flgstat = flagstat_metrics[sample]
        mkdup = markdups_metrics[sample]
        insrt = insert_metrics[sample]
        wgs = wgs_metrics[sample]
        outfile = os.path.join(tempdir, sample + "_metrics.csv")
        sample_outputs.append(outfile)

        collmet = CollectMetrics(wgs, insrt, flgstat,
                                 mkdup, outfile, sample)
        collmet.main()

    merged_metrics_csv = os.path.join(tempdir, 'merged_alignment_metrics.csv')
    csvutils.concatenate_csv(sample_outputs, merged_metrics_csv)

    hdfutils.convert_csv_to_hdf(
        merged_metrics_csv,
        merged_metrics,
        '/alignment/metrics/')

