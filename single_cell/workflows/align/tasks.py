from __future__ import division

import os
import logging

from single_cell.utils import bamutils
from single_cell.utils import csvutils
from single_cell.utils import helpers
from single_cell.utils import picardutils
from single_cell.utils.singlecell_copynumber_plot_utils import PlotMetrics
from single_cell.workflows.align.dtypes import dtypes

from .scripts import CollectMetrics
from .scripts import GenerateCNMatrix
from .scripts import SummaryMetrics


class LibraryContaminationError(Exception):
    pass


def add_contamination_status(
        infile, outfile,
        reference='grch37', ref_threshold=0.6, alt_threshold=0.2,
        strict_validation=True
):
    data = csvutils.read_csv_and_yaml(infile)

    data = data.set_index('cell_id', drop=False)

    fastqscreen_cols = [col for col in data.columns.values if col.startswith('fastqscreen_')]

    reference = "fastqscreen_{}".format(reference)
    if reference not in fastqscreen_cols:
        raise Exception("Could not find the fastq screen counts")

    alts = [col for col in fastqscreen_cols if not col == reference]

    data['is_contaminated'] = False

    perc_ref = data[reference] / data['total_reads']
    data.loc[perc_ref <= ref_threshold, 'is_contaminated'] = True

    for altcol in alts:
        perc_alt = data[altcol] / data['total_reads']
        data.loc[perc_alt > alt_threshold, 'is_contaminated'] = True

    col_type = dtypes()['metrics']['is_contaminated']
    data['is_contaminated'] = data['is_contaminated'].astype(col_type)

    csvutils.write_dataframe_to_csv_and_yaml(
        data, outfile, dtypes()['metrics'], write_header=True
    )

    # get cells that are contaminated and have enopugh human reads
    check_df = data.loc[data['is_contaminated'] == True]
    check_df['perc_ref'] = data[reference] / data['total_reads']
    check_df = check_df[check_df['perc_ref'] > ref_threshold]
    if strict_validation and (len(check_df) / len(data) > 0.2):
        logging.error("over 20% of cells are contaminated")


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
