import pypeliner
import pypeliner.managed as mgd
import pandas as pd
import glob
import re
import os
from single_cell.workflows.qc import tasks


def create_patient_workflow(pseudobulk_group, mafs, sample_all_snv_csvs,  
    mutationreport, merged_maf, high_impact_maf, merged_snvs, 
    merged_high_impact_snvs
):

    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1, }
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.transform(
        name='merge_mafs',
        func='single_cell.workflows.qc.tasks.merge_mafs',
        args=(
            mafs,
            mgd.OutputFile(merged_maf),
        ),
        kwargs = {"id_colname": True}
    )
    workflow.transform(
        name='filter_merged_maf',
        func='single_cell.workflows.qc.tasks.filter_maf_for_high_impact',
        args=(
            mgd.InputFile(merged_maf),
            mgd.OutputFile(high_impact_maf),
        ),
    )
    workflow.transform(
        name='merge_snvs',
        func='single_cell.workflows.qc.tasks.merge_snvs',
        args=(
            sample_all_snv_csvs,
            mgd.OutputFile(merged_snvs),
        ),
        kwargs = {"id_colname": True}
    )
    workflow.transform(
        name='filter_snvs',
        func='single_cell.workflows.qc.tasks.filter_snvs_for_high_impact',
        args=(
            mgd.InputFile(merged_snvs),
            mgd.OutputFile(merged_high_impact_snvs),
        )
    )

    workflow.transform(
        name='mutationreport',
        func='single_cell.workflows.qc.tasks.create_mutation_report',
        args=(
            pseudobulk_group,
            mgd.InputFile(merged_maf),
            mgd.InputFile(high_impact_maf),
            mgd.InputFile(merged_high_impact_snvs),
            mgd.OutputFile(mutationreport),
        ),
    )

    return workflow


def create_sample_level_plots(patient, cell_id, library_id, mappability_file, 
    strelka_file, museq_file, cosmic_status_file, snpeff_file, dbsnp_status_file, 
    trinuc_file, counts_file, destruct_breakpoint_annotation, destruct_breakpoint_counts, 
    lumpy_breakpoint_annotation, lumpy_breakpoint_evidence,
    haplotype_allele_data, annotation_metrics, hmmcopy_reads, hmmcopy_segs, 
    hmmcopy_metrics, alignment_metrics, gc_metrics, indel_file, reporthtml, maf, 
    snvs_all_csv, out_dir, config
):

    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1, }

    vep_refdir = config['vep_reference_dir']
    scp_qc_docker = config["docker"]

    prefix = os.path.join(out_dir,  patient, cell_id, library_id)

    mutations_per_cell_png =  os.path.join(prefix, "mutations_per_cell.png")
    summary_csv =  os.path.join(prefix, "summary.csv")
    snvs_high_impact_csv =  os.path.join(prefix, "snvs_high_impact.csv")
    trinuc_csv =  os.path.join(prefix, "trinuc.csv")
    snv_adjacent_distance_png =  os.path.join(prefix, "snv_adjacent_distance.png")
    snv_genome_count_png =  os.path.join(prefix,  "snv_genome_count.png")  
    snv_cell_counts_png =  os.path.join(prefix,  "snv_cell_counts.png")  
    snv_alt_counts_png =  os.path.join(prefix,  "snv_alt_counts.png")  

    rearranegementtype_distribution_destruct_unfiltered_png =  os.path.join(prefix,  
        "rearranegementtype_distribution_destruct_unfiltered.png"
    )  
    chromosome_types_destruct_unfiltered_png =  os.path.join(prefix, 
        "chromosome_types_destruct_unfiltered.png") 

    rearranegementtype_distribution_destruct_filtered_png =  os.path.join(prefix,  
        "rearranegementtype_distribution_destruct_filtered.png"
    )  
    chromosome_types_destruct_filtered_png =  os.path.join(prefix,  
        "chromosome_types_destruct_filtered.png") 

    rearranegementtype_distribution_lumpy_unfiltered_png =  os.path.join(prefix,  
        "rearranegementtype_distribution_lumpy_unfiltered.png"
    )  
    chromosome_types_lumpy_unfiltered_png =  os.path.join(prefix,  
        "chromosome_types_lumpy_unfiltered.png") 

    baf_plot_png =  os.path.join(prefix,  "BAFplot.png") 
    cn_plot_png =  os.path.join(prefix,  "CNplot.png") 
    datatype_summary_csv =  os.path.join(prefix,  "datatype_summary.csv")
    
    workflow = pypeliner.workflow.Workflow(ctx=ctx)
    
    workflow.transform(
        name='vcf2maf',
        func='single_cell.workflows.qc.tasks.vcf2maf',
        args=(
            mgd.InputFile(indel_file),
            mgd.OutputFile(maf),
            mgd.TempSpace('vcf2maf_temp'),
            vep_refdir,
        ),
        kwargs = (
            {"docker_image" : scp_qc_docker["vcf2maf"]}
        )
    )

    workflow.transform(
        name='qc_plots',
        func="single_cell.workflows.qc.scripts.single_cell_qc_plots.qc_plots",
        args=( 
            cell_id,
            mgd.InputFile(mappability_file),
            mgd.InputFile(strelka_file),
            mgd.InputFile(museq_file),
            mgd.InputFile(cosmic_status_file),
            mgd.InputFile(snpeff_file),
            mgd.InputFile(dbsnp_status_file),
            mgd.InputFile(trinuc_file),
            mgd.InputFile(counts_file),
            mgd.InputFile(destruct_breakpoint_annotation),
            mgd.InputFile(destruct_breakpoint_counts), 
            mgd.InputFile(lumpy_breakpoint_annotation),
            mgd.InputFile(lumpy_breakpoint_evidence),     
            mgd.InputFile(haplotype_allele_data),     
            mgd.InputFile(annotation_metrics),  
            mgd.InputFile(hmmcopy_reads),  
            mgd.InputFile(hmmcopy_segs),  
            mgd.InputFile(hmmcopy_metrics),  
            mgd.InputFile(alignment_metrics),  
            mgd.InputFile(gc_metrics),  
            library_id,
            prefix,
            mgd.OutputFile(mutations_per_cell_png), 
            mgd.OutputFile(summary_csv), 
            mgd.OutputFile(snvs_high_impact_csv), 
            mgd.OutputFile(snvs_all_csv), 
            mgd.OutputFile(trinuc_csv), 
            mgd.OutputFile(snv_adjacent_distance_png), 
            mgd.OutputFile(snv_genome_count_png), 
            mgd.OutputFile(snv_cell_counts_png), 
            mgd.OutputFile(snv_alt_counts_png), 
            mgd.OutputFile(rearranegementtype_distribution_destruct_unfiltered_png), 
            mgd.OutputFile(chromosome_types_destruct_unfiltered_png),
            mgd.OutputFile(rearranegementtype_distribution_destruct_filtered_png), 
            mgd.OutputFile(chromosome_types_destruct_filtered_png),
            mgd.OutputFile(rearranegementtype_distribution_lumpy_unfiltered_png), 
            mgd.OutputFile(chromosome_types_lumpy_unfiltered_png),  
            mgd.OutputFile(baf_plot_png), 
            mgd.OutputFile(cn_plot_png), 
            mgd.OutputFile(datatype_summary_csv),
        ),
    )

    workflow.transform(
        name='create_main_report',
        func="single_cell.workflows.qc.tasks.sample_level_report",
        args=(
            mgd.InputFile(mutations_per_cell_png), 
            mgd.InputFile(summary_csv), 
            mgd.InputFile(snvs_high_impact_csv), 
            mgd.InputFile(snvs_all_csv), 
            mgd.InputFile(trinuc_csv), 
            mgd.InputFile(snv_adjacent_distance_png), 
            mgd.InputFile(snv_genome_count_png), 
            mgd.InputFile(snv_cell_counts_png), 
            mgd.InputFile(snv_alt_counts_png), 
            mgd.InputFile(rearranegementtype_distribution_destruct_unfiltered_png), 
            mgd.InputFile(chromosome_types_destruct_unfiltered_png),
            mgd.InputFile(rearranegementtype_distribution_destruct_filtered_png), 
            mgd.InputFile(chromosome_types_destruct_filtered_png),
            mgd.InputFile(rearranegementtype_distribution_lumpy_unfiltered_png), 
            mgd.InputFile(chromosome_types_lumpy_unfiltered_png),  
            mgd.InputFile(baf_plot_png), 
            mgd.InputFile(cn_plot_png), 
            mgd.InputFile(datatype_summary_csv),
            mgd.InputFile(maf),
            mgd.OutputFile(reporthtml),
            out_dir,
            cell_id + "_" + library_id,
        ),
    )

    return workflow
