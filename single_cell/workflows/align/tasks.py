from __future__ import division

import os

from single_cell.utils import bamutils
from single_cell.utils import csvutils
from single_cell.utils import helpers
from single_cell.utils import picardutils
from single_cell.utils.singlecell_copynumber_plot_utils import PlotMetrics
from single_cell.workflows.align.dtypes import dtypes

from .scripts import CollectMetrics
from .scripts import GenerateCNMatrix
from .scripts import SummaryMetrics


def plot_metrics(metrics, output, plot_title, gc_matrix, gc_content):
    plot = PlotMetrics(
        metrics, output, plot_title, gcbias_matrix=gc_matrix,
        gc_content=gc_content, tablename='/alignment/metrics',
        gc_tablename='/alignment/gc_metrics')
    plot.plot_alignment_metrics()


def get_summary_metrics(infile, output):
    summ = SummaryMetrics(infile, output)
    summ.main()


def collect_gc(infiles, outfile, tempdir):
    helpers.makedirs(tempdir)

    tempouts = []
    for cell_id, infile in infiles.items():
        tempout = os.path.join(
            tempdir,
            "{}.parsed.csv.gz".format(cell_id))
        tempouts.append(tempout)
        gen_gc = GenerateCNMatrix(infile, tempout, ',',
                                  'NORMALIZED_COVERAGE', cell_id,
                                  'gcbias', dtypes()['gc'])
        gen_gc.main()

    csvutils.concatenate_csv(tempouts, outfile)


def collect_metrics(flagstat_metrics, markdups_metrics, insert_metrics,
                    wgs_metrics, tempdir, merged_metrics):
    helpers.makedirs(tempdir)
    sample_outputs = []

    for sample in flagstat_metrics.keys():
        flgstat = flagstat_metrics[sample]
        mkdup = markdups_metrics[sample]
        insrt = insert_metrics[sample]
        wgs = wgs_metrics[sample]
        outfile = os.path.join(tempdir, sample + "_metrics.csv.gz")
        sample_outputs.append(outfile)

        collmet = CollectMetrics(wgs, insrt, flgstat,
                                 mkdup, outfile, sample, dtypes()['metrics'])
        collmet.main()

    csvutils.concatenate_csv(sample_outputs, merged_metrics)


def picard_wgs_dup(
        input_bam, markdups_bam, markdups_metrics, tempdir,
        ref_genome, wgs_metrics, picard_wgs_params,
        picard_docker=None

):
    tempdir_markdups = os.path.join(tempdir, 'markdups')
    helpers.makedirs(tempdir_markdups)

    picardutils.bam_markdups(
        input_bam,
        markdups_bam,
        markdups_metrics,
        tempdir_markdups,
        docker_image=picard_docker
    )

    tempdir_wgs = os.path.join(tempdir, 'wgs')
    helpers.makedirs(tempdir_wgs)

    picardutils.bam_collect_wgs_metrics(
        input_bam,
        ref_genome,
        wgs_metrics,
        picard_wgs_params,
        tempdir_wgs,
        docker_image=picard_docker
    )


def picard_insert_gc_flagstat(
        input_bam, ref_genome, gc_metrics, gc_metrics_summary, gc_metrics_pdf,
        tempdir, flagstat_metrics, insert_metrics, insert_pdf, picard_docker=None,
        samtools_docker=None
):
    bamutils.bam_flagstat(
        input_bam,
        flagstat_metrics,
        docker_image=samtools_docker
    )

    gc_tempdir = os.path.join(tempdir, 'gc')
    helpers.makedirs(gc_tempdir)

    picardutils.bam_collect_gc_metrics(
        input_bam,
        ref_genome,
        gc_metrics,
        gc_metrics_summary,
        gc_metrics_pdf,
        gc_tempdir,
        docker_image=picard_docker
    )

    insert_tempdir = os.path.join(tempdir, 'insert')
    helpers.makedirs(insert_tempdir)
    picardutils.bam_collect_insert_metrics(
        input_bam,
        flagstat_metrics,
        insert_metrics,
        insert_pdf,
        insert_tempdir,
        docker_image=picard_docker
    )


def tar_align_data(infiles, tar_output, tempdir):
    helpers.makedirs(tempdir)

    for infile in infiles:
        for key, filepath in infile.items():
            temp_path = os.path.join(tempdir, '{}_{}'.format(key, os.path.basename(filepath)))
            helpers.shutil.copyfile(filepath, temp_path)

    helpers.make_tarfile(tar_output, tempdir)
