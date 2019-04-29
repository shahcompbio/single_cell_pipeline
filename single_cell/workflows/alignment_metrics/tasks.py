import os
from scripts import CollectMetrics
from scripts import GenerateCNMatrix
from scripts import SummaryMetrics

from single_cell.utils import picardutils
from single_cell.utils import bamutils
from single_cell.utils import helpers
from single_cell.utils import csvutils
from single_cell.utils import hdfutils

from single_cell.utils.singlecell_copynumber_plot_utils import PlotMetrics


def get_postprocess_metrics(infile, tempdir,
                    containers, markdups_metrics, flagstat_metrics):

    if not os.path.exists(tempdir):
        helpers.makedirs(tempdir)

    outfile = os.path.join(tempdir, 'markdps.bam')
    outfile_index = outfile + '.bai'

    picardutils.bam_markdups(infile, outfile, markdups_metrics, tempdir,
                             docker_image=containers['picard'])

    bamutils.bam_index(outfile, outfile_index, docker_image=containers['samtools'])
    bamutils.bam_flagstat(outfile, flagstat_metrics, docker_image=containers['samtools'])

def plot_metrics(metrics, output, plot_title, gc_matrix, gc_content):

    plot = PlotMetrics(
        metrics, output, plot_title, gcbias_matrix=gc_matrix,
        gc_content=gc_content, tablename='/alignment/metrics',
        gc_tablename='/alignment/gc_metrics')
    plot.plot_alignment_metrics()


def get_summary_metrics(infile, output):
    summ = SummaryMetrics(infile, output)
    summ.main()


def annotate_metrics(infile, sample_info, outfile, yamlfile=None):

    csvutils.annotate_metrics(infile, sample_info, outfile, yamlfile=yamlfile)


def merge_all_metrics(infiles, outfile):
    csvutils.merge_csv(infiles, outfile, "outer", "cell_id")


def bam_collect_wgs_metrics(
        bam_filename, ref_genome, metrics_filename, containers, picard_wgs_params, tempdir):

    picardutils.bam_collect_wgs_metrics(
        bam_filename,
        ref_genome,
        metrics_filename,
        picard_wgs_params,
        tempdir,
        docker_image=containers['picard'])


def bam_collect_gc_metrics(
        bam_filename, ref_genome, metrics_filename, summary_filename, chart_filename, tempdir, containers):
    picardutils.bam_collect_gc_metrics(
        bam_filename,
        ref_genome,
        metrics_filename,
        summary_filename,
        chart_filename,
        tempdir,
        docker_image=containers['picard'])



def bam_collect_insert_metrics(
        bam_filename, flagstat_metrics_filename, metrics_filename, histogram_filename, tempdir,containers):

    picardutils.bam_collect_insert_metrics(
        bam_filename,
        flagstat_metrics_filename,
        metrics_filename,
        histogram_filename,
        tempdir, docker_image=containers['picard'])


def collect_gc(infiles, outfile, tempdir, yamlfile=None):

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

    csvutils.concatenate_csv(tempouts, outfile)


def collect_metrics(flagstat_metrics, markdups_metrics, insert_metrics,
                    wgs_metrics, tempdir, merged_metrics, biobloom_count_metrics):

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
                                 mkdup, outfile, sample, biobloom_count_metrics[sample])
        collmet.main()

    csvutils.concatenate_csv(sample_outputs, merged_metrics)


