import pypeliner
import pypeliner.managed as mgd
import pandas as pd
import glob
import re
import os
from single_cell.workflows.qc import tasks

def create_pseudobulk_group_workflow(pseudobulk_group, mafs, sample_all_snv_csvs,  mutationreport, merged_maf, high_impact_maf, merged_snvs, merged_high_impact_snvs):

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


def create_sample_level_plots(cell_id, library_id, mappability_file, 
    strelka_file, museq_file, cosmic_status_file, snpeff_file, dbsnp_status_file, trinuc_file, 
    counts_file, breakpoint_annotation, breakpoint_counts, haplotype_allele_data, annotation_metrics,
    hmmcopy_reads, hmmcopy_segs, hmmcopy_metrics, alignment_metrics, gc_metrics, indel_file, reporthtml, maf, 
    snvs_all_csv, tmp_dir, out_dir, outpath
):

    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1, }

    tantalus = "/work/shah/tantalus/"
    prefix = os.path.join(tmp_dir,  outpath)
    outprefix = os.path.join(out_dir,  outpath)

    mutations_per_cell_png =  os.path.join(prefix, "mutations_per_cell.png")
    summary_csv =  os.path.join(prefix, "summary.csv")
    snvs_high_impact_csv =  os.path.join(prefix, "snvs_high_impact.csv")
    # snvs_all_csv =  os.path.join(prefix, "snvs_all.csv")
    trinuc_csv =  os.path.join(prefix, "trinuc.csv")
    snv_adjacent_distance_png =  os.path.join(prefix, "snv_adjacent_distance.png")
    snv_genome_count_png =  os.path.join(prefix,  "snv_genome_count.png")  
    snv_cell_counts_png =  os.path.join(prefix,  "snv_cell_counts.png")  
    snv_alt_counts_png =  os.path.join(prefix,  "snv_alt_counts.png")  
    rearranegementtype_distribution_png =  os.path.join(prefix,  "rearranegementtype_distribution.png")  
    chromosome_types_png =  os.path.join(prefix,  "chromosome_types.png") 
    BAFplot_png =  os.path.join(prefix,  "BAFplot.png") 
    CNplot_png =  os.path.join(prefix,  "CNplot.png") 
    datatype_summary_csv =  os.path.join(prefix,  "datatype_summary.csv")
    
    #hardcoded for now
    vepdata = "/work/shah/reference/vep/"
    genomeref = "/work/shah/reference/genomes/GRCh37-lite/GRCh37-lite.fa"

    workflow = pypeliner.workflow.Workflow(ctx=ctx)
    

    # if len(indelvcfs) > 1:
    #     for indelvcf, outputmaf in zip(indelvcfs, outputmafs):
    #         workflow.transform(
    #             name='vcf2maf',
    #             func='single_cell.workflows.qc.tasks.vcf2maf',
    #             args=(
    #                 mgd.InputFile(indelvcf),
    #                 mgd.TempOutputFile(outputmaf),
    #                 mgd.TempSpace('vcf2maf_temp'),
    #                 genomeref,
    #                 vepdata,
    #             ),
    #         )
    #     workflow.transform(
    #         name='merge_mafs',
    #         func='single_cell.workflows.qc.tasks.merge_mafs',
    #         args=(
    #             outputmafs,
    #             mgd.OutputFile(maf),
    #         ),
    #     )

    # else:

    workflow.transform(
        name='vcf2maf',
        func='single_cell.workflows.qc.tasks.vcf2maf',
        args=(
            mgd.InputFile(indel_file),
            mgd.OutputFile(maf),
            mgd.TempSpace('vcf2maf_temp'),
            genomeref,
            vepdata,
        ),
    )

    workflow.transform(
        name='scgenome_plots',
        func="single_cell.workflows.qc.scripts.scgenome-analysis.scgenome_analysis",
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
            mgd.InputFile(breakpoint_annotation),
            mgd.InputFile(breakpoint_counts),      
            mgd.InputFile(haplotype_allele_data),     
            mgd.InputFile(annotation_metrics),  
            mgd.InputFile(hmmcopy_reads),  
            mgd.InputFile(hmmcopy_segs),  
            mgd.InputFile(hmmcopy_metrics),  
            mgd.InputFile(alignment_metrics),  
            mgd.InputFile(gc_metrics),  
            library_id,
            tantalus,
            prefix,
            outprefix,
            mgd.OutputFile(mutations_per_cell_png), 
            mgd.OutputFile(summary_csv), 
            mgd.OutputFile(snvs_high_impact_csv), 
            mgd.OutputFile(snvs_all_csv), 
            mgd.OutputFile(trinuc_csv), 
            mgd.OutputFile(snv_adjacent_distance_png), 
            mgd.OutputFile(snv_genome_count_png), 
            mgd.OutputFile(snv_cell_counts_png), 
            mgd.OutputFile(snv_alt_counts_png), 
            mgd.OutputFile(rearranegementtype_distribution_png), 
            mgd.OutputFile(chromosome_types_png),
            mgd.OutputFile(BAFplot_png), 
            mgd.OutputFile(CNplot_png), 
            mgd.OutputFile(datatype_summary_csv),

        ),
    )

    mutations_per_cell_png =  os.path.join(outprefix, "mutations_per_cell.png")
    summary_csv =  os.path.join(outprefix, "summary.csv")
    snvs_high_impact_csv =  os.path.join(outprefix, "snvs_high_impact.csv")
    # snvs_all_csv =  os.path.join(outprefix, "snvs_all.csv")
    trinuc_csv =  os.path.join(outprefix, "trinuc.csv")
    snv_adjacent_distance_png =  os.path.join(outprefix, "snv_adjacent_distance.png")
    snv_genome_count_png =  os.path.join(outprefix,  "snv_genome_count.png")  
    snv_cell_counts_png =  os.path.join(outprefix,  "snv_cell_counts.png")  
    snv_alt_counts_png =  os.path.join(outprefix,  "snv_alt_counts.png")  
    rearranegementtype_distribution_png =  os.path.join(outprefix,  "rearranegementtype_distribution.png")  
    chromosome_types_png =  os.path.join(outprefix,  "chromosome_types.png") 
    BAFplot_png =  os.path.join(outprefix,  "BAFplot.png") 
    CNplot_png =  os.path.join(outprefix,  "CNplot.png") 
    datatype_summary_csv =  os.path.join(outprefix,  "datatype_summary.csv")

    workflow.transform(
        name='create_main_report',
        # ctx={'io': 1, 'mem': 8, 'disk': 100},
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
            mgd.InputFile(rearranegementtype_distribution_png), 
            mgd.InputFile(chromosome_types_png),
            mgd.InputFile(BAFplot_png), 
            mgd.InputFile(CNplot_png), 
            mgd.InputFile(datatype_summary_csv),
            mgd.InputFile(maf),
            mgd.OutputFile(reporthtml),
            out_dir,
            cell_id + "_" + library_id,


        ),
    )

    return workflow
