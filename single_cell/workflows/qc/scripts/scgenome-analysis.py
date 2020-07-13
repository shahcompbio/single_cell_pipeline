def scgenome_analysis(sample_id, ticket, snv_genotyping_ticket, library_id, dir, prefix, out_dir, mutations_per_cell, summary, 
                     snvs_high_impact, snvs_all, trinuc, snv_adjacent_distance, snv_genome_count, 
                     snv_cell_counts, snv_alt_counts, rearranegementtype_distribution, chromosome_types,
                     BAFplot, CNplot, datatype_summary):


# def scgenome_analysis(sample_id, ticket, snv_genotyping_ticket, library_id, dir, prefix, out_dir, CNplot):  
    import matplotlib
    matplotlib.use('Agg')
    from scgenome.snvdata import filter_snv_data
    from shutil import copyfile

    import argparse
    import scgenome.utils
    import csv

    import logging

    import os
    import sys
    import logging
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd

    import seaborn
    import numpy as np
    import pandas as pd
    import pylab
    import sklearn.preprocessing

    import scgenome
    import scgenome.db.qc
    import scgenome.cnplot
    import scgenome.cnfilter
    import scgenome.cnclones

    from scgenome.loaders.qc import load_qc_data
    import scgenome.loaders.snv
    
    #get rid of .tmp
    mutations_per_cell = mutations_per_cell[:-4]
    summary=summary[:-4]
    snvs_high_impact=snvs_high_impact[:-4]
    snvs_all=snvs_all[:-4]
    trinuc=trinuc[:-4]
    snv_adjacent_distance=snv_adjacent_distance[:-4] 
    snv_genome_count=snv_genome_count[:-4]
    snv_cell_counts=snv_cell_counts[:-4] 
    snv_alt_counts=snv_alt_counts[:-4] 
    rearranegementtype_distribution=rearranegementtype_distribution[:-4] 
    chromosome_types=chromosome_types[:-4] 
    BAFplot=BAFplot[:-4] 
    CNplot=CNplot[:-4]
    datatype_summary=datatype_summary[:-4]



    snv_results = scgenome.loaders.snv.load_snv_data(
        dir + ticket,
    )
        
    snv_data = snv_results["snv_data"]
    logging.info(dir + snv_genotyping_ticket)

    
    snv_results = scgenome.loaders.snv.load_snv_data(
        dir + snv_genotyping_ticket, positions=snv_data[['chrom', 'coord', 'ref', 'alt']].drop_duplicates(),
        filter_sample_id=sample_id, filter_library_id=library_id,
    )
    snv_count_data = snv_results["snv_count_data"]


    prefix = prefix + "/"

    filtered_data = filter_snv_data(
            snv_data,
            snv_count_data,
            num_cells_threshold=None,
            sum_alt_threshold=None,
            figures_prefix=prefix,
        )

    snv_data = filtered_data[0]
    snv_count_data = filtered_data[1]

    from scgenome.snvdata import run_bulk_snv_analysis
    allcells = snv_count_data.loc[:, ["cell_id"]].drop_duplicates()
    run_bulk_snv_analysis(snv_data, snv_count_data, allcells, results_prefix=prefix)

    #Plot number of mutations per cell
    percell = snv_count_data[snv_count_data.alt_counts > 0].groupby('cell_id').size().reset_index(name='Number of mutations')
    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_subplot(1, 1, 1)
    seaborn.distplot(percell["Number of mutations"])
    ax.set_title(f'Number of mutations per cell')
    fig.savefig(mutations_per_cell, bbox_inches='tight')

    #Write csv file with number of mutations and number of cells
    nmuts = len(snv_data[["chrom", "coord", "ref", "alt"]].drop_duplicates()["coord"])
    ncells = len(snv_count_data.cell_id.unique())
    data = {'ncells':[ncells], 'nmutations':[nmuts]}
    df = pd.DataFrame(data)
    df.to_csv(summary)

    # Get high impact mutations and write to file
    high_impact = (snv_data.query('effect_impact == "HIGH" | effect == "NON_SYNONYMOUS_CODING" | effect == "STOP_GAINED"')
        [[
            'chrom', 'coord', 'ref', 'alt',
            'gene_name', 'effect', 'effect_impact',
            'is_cosmic', 'max_museq_score', 'max_strelka_score',
            'amino_acid_change', 'num_cells',
        ]]
        .drop_duplicates())
    high_impact = high_impact.merge(
        snv_count_data.groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['alt_counts'].sum().reset_index())
    high_impact = high_impact.merge(
        snv_count_data.groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['ref_counts'].sum().reset_index())
    high_impact = high_impact.merge(
        snv_count_data.groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['total_counts'].sum().reset_index())

    high_impact = high_impact.sort_values(by=['is_cosmic', "gene_name", 'num_cells'])
    high_impact.to_csv(snvs_high_impact)

    snv_data = snv_data.merge(
        snv_count_data.groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['alt_counts'].sum().reset_index())
    snv_data = snv_data.merge(
        snv_count_data.groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['ref_counts'].sum().reset_index())
    snv_data = snv_data.merge(
        snv_count_data.groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['total_counts'].sum().reset_index())

    snv_data.to_csv(snvs_all)

    # Add extra columns for number of cells class
    snv_data = snv_data[snv_data['num_cells'] > 0]
    snv_data['num_cells_class'] = '1'
    snv_data.loc[snv_data['num_cells'] == 1, 'num_cells_class'] = '1'
    snv_data.loc[snv_data['num_cells'] > 1, 'num_cells_class'] = '2-5'
    snv_data.loc[snv_data['num_cells'] > 5, 'num_cells_class'] = '6-20'
    snv_data.loc[snv_data['num_cells'] > 20, 'num_cells_class'] = '>20'

    import wgs_analysis.snvs.mutsig
    import wgs_analysis.plots.snv
    import wgs_analysis.annotation.position

    #Sum number of mutation per 96 channel
    sigs, sig_prob = wgs_analysis.snvs.mutsig.load_signature_probabilities()
    snv_data_sig = snv_data[snv_data['tri_nucleotide_context'].notnull()]
    snv_data_sig = snv_data_sig.merge(sigs)
    tri_nuc_table = (snv_data_sig.groupby(['num_cells_class', 'tri_nuc_idx'])
                    .size().reset_index().merge(sigs[0:95])
                    .rename(columns={0: "n"}))
    tri_nuc_table.to_csv(trinuc)

    ##### Make PNG plots

    # Annotate adjacent distance
    snv_data = wgs_analysis.annotation.position.annotate_adjacent_distance(snv_data)

    # Plot adjacent distance of SNVs
    fig = plt.figure(figsize=(10, 3))
    wgs_analysis.plots.snv.snv_adjacent_distance_plot(plt.gca(), snv_data)
    fig.savefig(snv_adjacent_distance, bbox_inches='tight')

    # Plot snv count as a histogram across the genome
    fig = plt.figure(figsize=(10, 3))
    wgs_analysis.plots.snv.snv_count_plot(plt.gca(), snv_data)
    fig.savefig(snv_genome_count, bbox_inches='tight')

    # Calculate cell counts
    cell_counts = (
        snv_count_data
        .query('alt_counts > 0')
        .drop_duplicates(['chrom', 'coord', 'ref', 'alt', 'cell_id'])
        .groupby(['chrom', 'coord', 'ref', 'alt'], observed=True).size().rename('num_cells')
        .reset_index())

    fig = plt.figure(figsize=(4, 4))
    cell_counts['num_cells'].astype(float).hist(bins=50)
    fig.savefig(snv_cell_counts, bbox_inches='tight')

    snv_data = snv_data.merge(cell_counts, how='left')
    assert not snv_data['num_cells'].isnull().any()

    # Calculate total alt counts for each SNV
    sum_alt_counts = (
        snv_count_data
        .groupby(['chrom', 'coord', 'ref', 'alt'], observed=True)['alt_counts']
        .sum().rename('sum_alt_counts')
        .reset_index())

    fig = plt.figure(figsize=(4, 4))
    sum_alt_counts['sum_alt_counts'].astype(float).hist(bins=50)
    fig.savefig(snv_alt_counts, bbox_inches='tight')


    # ########################################
    # # Analyse breakpoint data
    # ########################################

    import scgenome.loaders.breakpoint

    breakpoint_results = scgenome.loaders.breakpoint.load_breakpoint_data(
        os.path.join(dir, ticket)
    )

    import scgenome.breakpointdata

    breakpoint_data = breakpoint_results['breakpoint_data']
    breakpoint_data["library_id"] = library_id
    breakpoint_data["sample_id"] = sample_id
    breakpoint_count_data = breakpoint_results['breakpoint_count_data']
    breakpoint_count_data["library_id"] = library_id
    breakpoint_count_data["sample_id"] = sample_id

    breakpoint_data, breakpoint_count_data = scgenome.breakpointdata.annotate_breakpoint_data(
        breakpoint_data,
        breakpoint_count_data,
    )

    breakpoint_data, breakpoint_count_data = scgenome.breakpointdata.filter_breakpoint_data(
        breakpoint_data,
        breakpoint_count_data,
    )


    # Plot distribution of breakpoint types

    features = [
        ('log_num_reads', 'log discordant read counts'),
        ('log_num_split', 'log split read counts'),
        ('template_length_min', 'prediction sequence length'),
        ('homology', 'sequence homology'),
        ('log_num_cells', 'log number of cells'),
    ]
    print("wqokncon")
    if len(breakpoint_data.sample_id) == 0:
        fig = plt.figure(figsize=(8, 12))
        fig.savefig(rearranegementtype_distribution, bbox_inches='tight')
        #fig.savefig(prefix + 'breakpoint_adjacent_distance.png', bbox_inches='tight')
        fig.savefig(chromosome_types, bbox_inches='tight')
    elif (len(breakpoint_data.sample_id) < 20) & (len(breakpoint_data.sample_id) > 0):
        order = breakpoint_data['library_id'].unique()
        hue_order = breakpoint_data['rearrangement_type'].unique()

        fig = plt.figure(figsize=(8, 12))
        ax = fig.add_subplot(len(features) + 1, 1, 1)
        plot_data = (
            breakpoint_data
            .groupby(['rearrangement_type'])
            .size().rename('count').reset_index())
        seaborn.barplot(
            ax=ax, data=plot_data,
            x='rearrangement_type', y='count',
            hue_order=hue_order)
        ax.set_title(f'Counts by rearrangement type')
        #ax.set_xticklabels([])
        ax.set_xlabel('', visible=False)
        ax.set_ylabel('Counts', visible=False)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        fig.savefig(rearranegementtype_distribution, bbox_inches='tight')

        fig = plt.figure(figsize=(8, 12))
        fig.savefig(chromosome_types, bbox_inches='tight')
    else:
        order = breakpoint_data['library_id'].unique()
        hue_order = breakpoint_data['rearrangement_type'].unique()

        fig = plt.figure(figsize=(8, 12))
        ax = fig.add_subplot(len(features) + 1, 1, 1)
        plot_data = (
            breakpoint_data
            .groupby(['rearrangement_type'])
            .size().rename('count').reset_index())
        seaborn.barplot(
            ax=ax, data=plot_data,
            x='rearrangement_type', y='count',
            hue_order=hue_order)
        ax.set_title(f'Counts by rearrangement type')
        #ax.set_xticklabels([])
        ax.set_xlabel('', visible=False)
        ax.set_ylabel('Counts', visible=False)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        fig.savefig(rearranegementtype_distribution, bbox_inches='tight')


        # Plot rearrangement type distribution across the genome per library
        fig = plt.figure(figsize=(10, 3))
        breakends = wgs_analysis.plots.rearrangement.create_breakends(
            breakpoint_data, data_cols=['rearrangement_type'])
        ax = fig.add_subplot(1, 1, 1)
        wgs_analysis.plots.rearrangement.chromosome_type_plot(
            ax, breakends)
        ax.set_title(f'Chromosome types')
        plt.tight_layout()
        fig.savefig(chromosome_types, bbox_inches='tight')


    # ########################################
    # # Analyse SNP data
    # # ########################################

    import scgenome.loaders.allele
    allele_results = scgenome.loaders.allele.load_haplotype_allele_data(
        dir + ticket,
    )
    allele_data = allele_results['allele_counts']
    index_cols = [
        'chromosome',
        'start',
        'end',
        'hap_label',
    ]
    allele_data = allele_data.set_index(index_cols + ['cell_id', 'allele_id'])['readcount'].unstack(fill_value=0)
    allele_data.rename(columns={0: 'allele_1', 1: 'allele_2'}, inplace=True)
    allele_data.reset_index(inplace=True)
    allele_data['total_counts'] = allele_data['allele_1'] + allele_data['allele_2']
    allele_data = allele_data.astype({"chromosome": "str", "cell_id":"str"})

    merged_allele_data = (
        allele_data.groupby(['chromosome', 'hap_label', 'end', 'start'])[['allele_1', 'allele_2', 'total_counts']]
        .aggregate({'allele_1': sum, 'allele_2': sum, 'total_counts': sum}).reset_index())

    merged_allele_data['baf'] = np.minimum(merged_allele_data['allele_1'], merged_allele_data['allele_2']) / merged_allele_data['total_counts'].astype(float)
    import scgenome.snpdata
    merged_allele_data = merged_allele_data.rename(columns={
            'chromosome': 'chr'})
    fig = plt.figure(figsize=(20, 4))
    ax = fig.add_subplot(1, 1, 1)
    scgenome.snpdata.plot_vaf_profile(
        ax, merged_allele_data, 'baf',
        size_field_name='total_counts',
        size_scale=5000.,
    )
    fig.savefig( BAFplot, bbox_inches='tight')


    # ########################################
    # # Plot CN profile
    # ########################################

    import scgenome.db.qc
    import scgenome.cnplot
    import scgenome.cnfilter
    import scgenome.cnclones

    results_tables = scgenome.db.qc.get_qc_data(
        [ticket],
        dir,
        sample_ids=[sample_id],
        do_caching=False,
    )

    cn_data, metrics_data = (
        results_tables['hmmcopy_reads'],
        results_tables['annotation_metrics'],
    )

    metrics_data = scgenome.cnfilter.calculate_filter_metrics(
        metrics_data,
        cn_data,
    )

    filtered_cells = metrics_data.loc[
        metrics_data['filter_is_s_phase'] &
        metrics_data['filter_quality'],
        ['cell_id']]

    logging.info('filtering {} of {} cells'.format(
        len(filtered_cells.index), len(metrics_data.index)))

    cn_data_filt = cn_data.merge(filtered_cells[['cell_id']].drop_duplicates())

    import statistics
    if len(cn_data_filt.chr) == 0:
        fig = plt.figure(figsize=(8, 12))
        fig.savefig( CNplot, bbox_inches='tight')
    else:
        cn_data_summary = (
            cn_data_filt.groupby(['chr', 'start', 'end'])[['copy', 'state']]
            .aggregate({'copy': statistics.median, 'state': statistics.median}).reset_index()).dropna()
        cn_data_summary['state'] = cn_data_summary['state'].astype(float).round().astype(int)
        import scgenome.cnplot
        fig = plt.figure(figsize=(20, 4))
        ax = fig.add_subplot(1, 1, 1)
        scgenome.cnplot.plot_cell_cn_profile(
            ax, cn_data_summary, 'copy', 'state'
        )
        fig.savefig(CNplot, bbox_inches='tight')

    #######################################
    #number of cells with each data type
    #######################################
    d = {'datatype': ['snv', 'cnv', 'haplotypes'], 'ncells': [len(snv_count_data.cell_id.unique()), len(metrics_data["cell_id"].unique()), len(allele_data["cell_id"].unique())]}
    df = pd.DataFrame(d)
    df.to_csv(datatype_summary)


    plots = [os.path.join(prefix, f) for f in os.listdir(prefix)]

    pseudobulk_group = sample_id + "_" + library_id
    os.makedirs(os.path.join(out_dir, pseudobulk_group))
    for plot in plots:
        # os.rename(plot, os.path.join(out_dir, pseudobulk_group, os.path.basename(plot)))
        copyfile(plot, os.path.join(out_dir, pseudobulk_group, os.path.basename(plot) + ".tmp"))
        copyfile(plot,plot+ ".tmp")
        copyfile(plot, os.path.join(out_dir, pseudobulk_group, os.path.basename(plot)))


    return

