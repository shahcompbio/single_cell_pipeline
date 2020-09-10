import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn
import numpy as np
import pandas as pd
import statistics
import wgs_analysis.plots.rearrangement
import wgs_analysis.snvs.mutsig
import wgs_analysis.plots.snv
import wgs_analysis.annotation.position
from scgenome.loaders.qc import load_qc_data
from scgenome.snvdata import run_bulk_snv_analysis
import scgenome
import scgenome.db.qc_from_files
import scgenome.cnplot
import scgenome.cnfilter
import scgenome.cnclones
import scgenome.loaders.allele
import scgenome.snpdata
import scgenome.loaders.snv
import scgenome.loaders.breakpoint
import scgenome.breakpointdata
import scgenome.utils
from scgenome.snvdata import filter_snv_data


def load_snv_data(
        sample_id, library_id, prefix, mappability_file, strelka_file,
        museq_file, cosmic_status_file, snpeff_file, dbsnp_status_file,
        trinuc_file, counts_file
):
    snv_results = scgenome.loaders.snv.load_snv_data_from_files(
        [mappability_file],
        [strelka_file],
        [museq_file],
        [cosmic_status_file],
        [snpeff_file],
        [dbsnp_status_file],
        [trinuc_file],
        [counts_file],
        snv_annotation=True
    )

    snv_data = snv_results["snv_data"]

    snv_results = scgenome.loaders.snv.load_snv_data_from_files(
        [mappability_file],
        [strelka_file],
        [museq_file],
        [cosmic_status_file],
        [snpeff_file],
        [dbsnp_status_file],
        [trinuc_file],
        [counts_file],
        snv_counts=True,
        positions=snv_data[['chrom', 'coord', 'ref', 'alt']].drop_duplicates(),
        filter_sample_id=sample_id,
        filter_library_id=library_id
    )

    snv_count_data = snv_results["snv_count_data"]

    filtered_data = filter_snv_data(
        snv_data,
        snv_count_data,
        num_cells_threshold=None,
        sum_alt_threshold=None,
        figures_prefix=prefix + "/",
    )

    return filtered_data[0], filtered_data[1]


def plot_mutations_per_cell(snv_data, snv_count_data, mutations_per_cell, prefix):
    allcells = snv_count_data.loc[:, ["cell_id"]].drop_duplicates()
    run_bulk_snv_analysis(snv_data, snv_count_data, allcells, results_prefix=prefix)

    # Plot number of mutations per cell
    percell = snv_count_data[snv_count_data.alt_counts > 0].groupby('cell_id').size()
    percell = percell.reset_index(name='Number of mutations')

    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_subplot(1, 1, 1)
    seaborn.distplot(percell["Number of mutations"])
    ax.set_title('Number of mutations per cell')
    fig.savefig(mutations_per_cell, bbox_inches='tight', format="png")


def write_summary_csv(snv_data, snv_count_data, summary):
    # Write csv file with number of mutations and number of cells
    nmuts = len(snv_data[["chrom", "coord", "ref", "alt"]].drop_duplicates()["coord"])
    ncells = len(snv_count_data.cell_id.unique())
    data = {'ncells': [ncells], 'nmutations': [nmuts]}
    df = pd.DataFrame(data)
    df.to_csv(summary)


def write_high_impact_snvs(snv_data, snv_count_data, snvs_high_impact, snvs_all):
    # Get high impact mutations and write to file
    filter_query = 'effect_impact == "HIGH"| effect == "NON_SYNONYMOUS_CODING" | effect == "STOP_GAINED"'
    relevant_cols = [
        'chrom', 'coord', 'ref', 'alt',
        'gene_name', 'effect', 'effect_impact',
        'is_cosmic', 'max_museq_score', 'max_strelka_score',
        'amino_acid_change', 'num_cells',
    ]

    high_impact = snv_data.query(filter_query)
    high_impact = high_impact[relevant_cols].drop_duplicates()

    snv_count_data = snv_count_data.groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)

    high_impact = high_impact.merge(snv_count_data['alt_counts'].sum().reset_index())
    high_impact = high_impact.merge(snv_count_data['ref_counts'].sum().reset_index())
    high_impact = high_impact.merge(snv_count_data['total_counts'].sum().reset_index())
    high_impact = high_impact.sort_values(by=['is_cosmic', "gene_name", 'num_cells'])
    high_impact.to_csv(snvs_high_impact)

    snv_data = snv_data.merge(snv_count_data['alt_counts'].sum().reset_index())
    snv_data = snv_data.merge(snv_count_data['ref_counts'].sum().reset_index())
    snv_data = snv_data.merge(snv_count_data['total_counts'].sum().reset_index())

    snv_data.to_csv(snvs_all)


def write_trinucleotide_context(snv_data, trinuc):
    # Sum number of mutation per 96 channel
    sigs, sig_prob = wgs_analysis.snvs.mutsig.load_signature_probabilities()
    snv_data_sig = snv_data[snv_data['tri_nucleotide_context'].notnull()]
    snv_data_sig = snv_data_sig.merge(sigs)

    tri_nuc_table = snv_data_sig.groupby(['num_cells_class', 'tri_nuc_idx'])
    tri_nuc_table = tri_nuc_table.size().reset_index().merge(sigs[0:95])
    tri_nuc_table = tri_nuc_table.rename(columns={0: "n"})

    tri_nuc_table.to_csv(trinuc)


def plot_snv_adjacent_distance(snv_data, snv_adjacent_distance):
    # Plot adjacent distance of SNVs
    fig = plt.figure(figsize=(10, 3))
    wgs_analysis.plots.snv.snv_adjacent_distance_plot(plt.gca(), snv_data)
    fig.savefig(snv_adjacent_distance, bbox_inches='tight', format="png")


def plot_snv_genome_count(snv_data, snv_genome_count):
    # Plot snv count as a histogram across the genome
    fig = plt.figure(figsize=(10, 3))
    wgs_analysis.plots.snv.snv_count_plot(plt.gca(), snv_data)
    fig.savefig(snv_genome_count, bbox_inches='tight', format="png")


def plot_snv_cell_counts(cell_counts, snv_cell_counts):
    fig = plt.figure(figsize=(4, 4))
    cell_counts['num_cells'].astype(float).hist(bins=50)
    fig.savefig(snv_cell_counts, bbox_inches='tight', format="png")


def plot_snv_alt_counts(snv_count_data, snv_alt_counts):
    # Calculate total alt counts for each SNV
    snv_count_data = snv_count_data.groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)
    sum_alt_counts = snv_count_data['alt_counts'].sum().rename('sum_alt_counts').reset_index()

    fig = plt.figure(figsize=(4, 4))
    sum_alt_counts['sum_alt_counts'].astype(float).hist(bins=50)
    fig.savefig(snv_alt_counts, bbox_inches='tight', format="png")


def load_breakpoint_data(
        sample_id, library_id, breakpoint_annotation, breakpoint_count,
        lumpy=False, filter_data=False
):
    breakpoint_results = scgenome.loaders.breakpoint.load_breakpoint_data_from_files(
        [breakpoint_annotation],
        [breakpoint_count],
        lumpy=lumpy
    )

    breakpoint_data = breakpoint_results['breakpoint_data']
    breakpoint_data["library_id"] = library_id
    breakpoint_data["sample_id"] = sample_id
    breakpoint_count_data = breakpoint_results['breakpoint_count_data']
    breakpoint_count_data["library_id"] = library_id
    breakpoint_count_data["sample_id"] = sample_id

    breakpoint_data, breakpoint_count_data = scgenome.breakpointdata.annotate_breakpoint_data(
        breakpoint_data,
        breakpoint_count_data,
        is_lumpy=lumpy
    )

    if filter_data:
        breakpoint_data, breakpoint_count_data = scgenome.breakpointdata.filter_breakpoint_data(
            breakpoint_data, breakpoint_count_data,
        )

    return breakpoint_data


def plot_breakpoint_distribution(
        breakpoint_data, rearranegementtype_distribution, chromosome_types, is_lumpy=False
):
    features = [
        ('log_num_reads', 'log discordant read counts'),
        ('log_num_split', 'log split read counts'),
        ('template_length_min', 'prediction sequence length'),
        ('homology', 'sequence homology'),
        ('log_num_cells', 'log number of cells'),
    ]
    type_col = "rearrangement_type"

    if is_lumpy:
        breakpoint_data = breakpoint_data.rename(
            columns={
                'breakpoint_id': "prediction_id",
                'chrom1': "chromosome_1", 'strand1': "strand_1",
                'start1': "position_1",
                'chrom2': "chromosome_2", 'strand2': "strand_2",
                'start2': "position_2", "type": "rearrangement_type"
            }
        )

    hue_order = breakpoint_data[type_col].unique()

    if len(breakpoint_data.sample_id) == 0:
        fig = plt.figure(figsize=(8, 12))
        fig.savefig(rearranegementtype_distribution, bbox_inches='tight', format="png")
        fig.savefig(chromosome_types, bbox_inches='tight', format="png")
    elif (len(breakpoint_data.sample_id) < 20) and (len(breakpoint_data.sample_id) > 0):

        fig = plt.figure(figsize=(8, 12))
        ax = fig.add_subplot(len(features) + 1, 1, 1)
        plot_data = (
            breakpoint_data.groupby([type_col]).size().rename('count').reset_index()
        )
        seaborn.barplot(
            ax=ax, data=plot_data, x=type_col, y='count', hue_order=hue_order
        )

        ax.set_title(f'Counts by rearrangement type')
        ax.set_xlabel('', visible=False)
        ax.set_ylabel('Counts', visible=False)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        fig.savefig(rearranegementtype_distribution, bbox_inches='tight', format="png")
        fig = plt.figure(figsize=(8, 12))
        fig.savefig(chromosome_types, bbox_inches='tight', format="png")
    else:
        fig = plt.figure(figsize=(8, 12))
        ax = fig.add_subplot(len(features) + 1, 1, 1)
        plot_data = breakpoint_data.groupby([type_col]).size().rename('count').reset_index()

        seaborn.barplot(
            ax=ax, data=plot_data, x=type_col, y='count', hue_order=hue_order
        )

        ax.set_title('Counts by rearrangement type')
        ax.set_xlabel('', visible=False)
        ax.set_ylabel('Counts', visible=False)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        fig.savefig(rearranegementtype_distribution, bbox_inches='tight', format="png")

        # Plot rearrangement type distribution across the genome per library
        fig = plt.figure(figsize=(10, 3))

        breakends = wgs_analysis.plots.rearrangement.create_breakends(
            breakpoint_data, data_cols=[type_col]
        )
        ax = fig.add_subplot(1, 1, 1)
        wgs_analysis.plots.rearrangement.chromosome_type_plot(
            ax, breakends, rearrangement_types=hue_order
        )
        ax.set_title(f'Chromosome types')
        plt.tight_layout()
        fig.savefig(chromosome_types, bbox_inches='tight', format="png")


def load_allele_data(haplotype_allele_data):
    allele_results = scgenome.loaders.allele.load_haplotype_allele_data_from_file([haplotype_allele_data])

    allele_data = allele_results['allele_counts']
    index_cols = [
        'chromosome', 'start', 'end', 'hap_label',
    ]
    allele_data = allele_data.set_index(index_cols + ['cell_id', 'allele_id'])['readcount'].unstack(fill_value=0)
    allele_data.rename(columns={0: 'allele_1', 1: 'allele_2'}, inplace=True)
    allele_data.reset_index(inplace=True)
    allele_data['total_counts'] = allele_data['allele_1'] + allele_data['allele_2']
    allele_data = allele_data.astype({"chromosome": "str", "cell_id": "str"})

    return allele_data


def plotbaf(allele_data, baf_plot):
    merged_allele_data = allele_data.groupby(['chromosome', 'hap_label', 'end', 'start'])
    merged_allele_data = merged_allele_data[['allele_1', 'allele_2', 'total_counts']]
    merged_allele_data = merged_allele_data.aggregate(
        {'allele_1': sum, 'allele_2': sum, 'total_counts': sum}
    )
    merged_allele_data = merged_allele_data.reset_index()

    baf_numerator = np.minimum(merged_allele_data['allele_1'], merged_allele_data['allele_2'])
    merged_allele_data['baf'] = baf_numerator / merged_allele_data['total_counts'].astype(float)

    merged_allele_data = merged_allele_data.rename(
        columns={'chromosome': 'chr'}
    )

    fig = plt.figure(figsize=(20, 4))
    ax = fig.add_subplot(1, 1, 1)

    scgenome.snpdata.plot_vaf_profile(
        ax, merged_allele_data, 'baf',
        size_field_name='total_counts',
        size_scale=5000.,
    )
    fig.savefig(baf_plot, bbox_inches='tight', format="png")


def load_qc_data(
        sample_id, annotation_metrics, hmmcopy_reads, hmmcopy_segs,
        hmmcopy_metrics, alignment_metrics, gc_metrics
):
    results_tables_new = scgenome.db.qc_from_files.get_qc_data_from_filenames(
        [annotation_metrics],
        [hmmcopy_reads],
        [hmmcopy_segs],
        [hmmcopy_metrics], [alignment_metrics], [gc_metrics],
        sample_ids=[sample_id], additional_hmmcopy_reads_cols=None
    )

    cn_data = results_tables_new['hmmcopy_reads']
    metrics_data = results_tables_new['annotation_metrics']

    metrics_data = scgenome.cnfilter.calculate_filter_metrics(
        metrics_data,
        cn_data,
    )

    filtered_cells = metrics_data.loc[
        metrics_data['filter_is_s_phase'] &
        metrics_data['filter_quality'],
        ['cell_id']]

    cn_data_filt = cn_data.merge(filtered_cells[['cell_id']].drop_duplicates())

    return metrics_data, cn_data_filt


def plot_cn(cn_data_filt, cn_plot):
    if len(cn_data_filt.chr) == 0:
        fig = plt.figure(figsize=(8, 12))
        fig.savefig(cn_plot, bbox_inches='tight', format="png")
    else:
        cn_data_summary = cn_data_filt.groupby(['chr', 'start', 'end'])
        cn_data_summary = cn_data_summary[['copy', 'state']]
        cn_data_summary = cn_data_summary.aggregate(
            {'copy': statistics.median, 'state': statistics.median}
        )
        cn_data_summary = cn_data_summary.reset_index().dropna()

        cn_data_summary['state'] = cn_data_summary['state'].astype(float).round().astype(int)

        fig = plt.figure(figsize=(20, 4))
        ax = fig.add_subplot(1, 1, 1)
        scgenome.cnplot.plot_cell_cn_profile(
            ax, cn_data_summary, 'copy', 'state'
        )
        fig.savefig(cn_plot, bbox_inches='tight', format="png")


def qc_plots(
        sample_id, mappability_file, strelka_file, museq_file,
        cosmic_status_file, snpeff_file, dbsnp_status_file, trinuc_file,
        counts_file, destruct_breakpoint_annotation, destruct_breakpoint_count,
        lumpy_breakpoint_annotation, lumpy_breakpoint_evidence, haplotype_allele_data,
        annotation_metrics, hmmcopy_reads, hmmcopy_segs, hmmcopy_metrics, alignment_metrics,
        gc_metrics, library_id, prefix, mutations_per_cell, summary, snvs_high_impact, snvs_all,
        trinuc, snv_adjacent_distance, snv_genome_count, snv_cell_counts, snv_alt_counts,
        rearranegementtype_distribution_destruct_unfiltered, chromosome_types_destruct_unfiltered,
        rearranegementtype_distribution_destruct_filtered, chromosome_types_destruct_filtered,
        rearranegementtype_distribution_lumpy_unfiltered, chromosome_types_lumpy_unfiltered,
        baf_plot, cn_plot, datatype_summary
):
    snv_data, snv_count_data = load_snv_data(
        sample_id, library_id, prefix, mappability_file, strelka_file, museq_file,
        cosmic_status_file, snpeff_file, dbsnp_status_file, trinuc_file, counts_file
    )

    plot_mutations_per_cell(snv_data, snv_count_data, mutations_per_cell, prefix)

    write_summary_csv(snv_data, snv_count_data, summary)

    write_high_impact_snvs(snv_data, snv_count_data, snvs_high_impact, snvs_all)

    # Add extra columns for number of cells class
    snv_data = snv_data[snv_data['num_cells'] > 0]
    snv_data['num_cells_class'] = '1'
    snv_data.loc[snv_data['num_cells'] == 1, 'num_cells_class'] = '1'
    snv_data.loc[snv_data['num_cells'] > 1, 'num_cells_class'] = '2-5'
    snv_data.loc[snv_data['num_cells'] > 5, 'num_cells_class'] = '6-20'
    snv_data.loc[snv_data['num_cells'] > 20, 'num_cells_class'] = '>20'

    write_trinucleotide_context(snv_data, trinuc)

    # Annotate adjacent distance
    snv_data = wgs_analysis.annotation.position.annotate_adjacent_distance(snv_data)

    plot_snv_adjacent_distance(snv_data, snv_adjacent_distance)

    plot_snv_genome_count(snv_data, snv_genome_count)

    cell_counts = snv_count_data.query('alt_counts > 0')
    cell_counts = cell_counts.drop_duplicates(['chrom', 'coord', 'ref', 'alt', 'cell_id'])
    cell_counts = cell_counts.groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)
    cell_counts = cell_counts.size().rename('num_cells').reset_index()

    plot_snv_cell_counts(cell_counts, snv_cell_counts)

    snv_data = snv_data.merge(cell_counts, how='left')
    assert not snv_data['num_cells'].isnull().any()

    plot_snv_alt_counts(snv_count_data, snv_alt_counts)

    # Analyse breakpoint data

    destruct_breakpoint_data_unfiltered = load_breakpoint_data(
        sample_id, library_id,
        destruct_breakpoint_annotation,
        destruct_breakpoint_count)

    plot_breakpoint_distribution(
        destruct_breakpoint_data_unfiltered,
        rearranegementtype_distribution_destruct_unfiltered,
        chromosome_types_destruct_unfiltered)

    destruct_breakpoint_data_filtered = load_breakpoint_data(
        sample_id, library_id,
        destruct_breakpoint_annotation, destruct_breakpoint_count,
        filter_data=True)

    plot_breakpoint_distribution(
        destruct_breakpoint_data_filtered,
        rearranegementtype_distribution_destruct_filtered,
        chromosome_types_destruct_filtered)

    lumpy_breakpoint_data_unfiltered = load_breakpoint_data(
        sample_id, library_id,
        lumpy_breakpoint_annotation, lumpy_breakpoint_evidence,
        lumpy=True)

    plot_breakpoint_distribution(
        lumpy_breakpoint_data_unfiltered,
        rearranegementtype_distribution_lumpy_unfiltered,
        chromosome_types_lumpy_unfiltered, is_lumpy=True)

    # Analyse SNP data
    allele_data = load_allele_data(haplotype_allele_data)
    plotbaf(allele_data, baf_plot)

    # Plot CN profile
    metrics_data, cn_data_filt = load_qc_data(
        sample_id, annotation_metrics, hmmcopy_reads,
        hmmcopy_segs, hmmcopy_metrics, alignment_metrics, gc_metrics)

    plot_cn(cn_data_filt, cn_plot)

    # number of cells with each data type
    df = {
        'datatype': ['snv', 'cnv', 'haplotypes'],
        'ncells': [
            len(snv_count_data.cell_id.unique()),
            len(metrics_data["cell_id"].unique()),
            len(allele_data["cell_id"].unique())
        ]
    }
    df = pd.DataFrame(df)
    df.to_csv(datatype_summary)
