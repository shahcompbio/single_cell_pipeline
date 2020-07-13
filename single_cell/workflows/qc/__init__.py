import pypeliner
import pypeliner.managed as mgd
import pandas as pd
import glob
import re
import os

def create_qc_workflow(cell_id, jira_id, snv_jira_id, library_id, tmp_dir, out_dir, reporthtml, mutationreporthtml, maf):
    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1, }

    tantalus = "/work/shah/tantalus/"
    prefix = os.path.join(tmp_dir,  cell_id + "_" + library_id)

    mutations_per_cell_png =  os.path.join(prefix, "mutations_per_cell.png")
    summary_csv =  os.path.join(prefix, "summary.csv")
    snvs_high_impact_csv =  os.path.join(prefix, "snvs_high_impact.csv")
    snvs_all_csv =  os.path.join(prefix, "snvs_all.csv")
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


    workflow = pypeliner.workflow.Workflow(ctx=ctx)
    print(mutations_per_cell_png)
    workflow.transform(
        name='scgenome_plots',
        # ctx={'io': 1, 'mem': 8, 'disk': 100},
        func="single_cell.workflows.qc.scripts.scgenome-analysis.scgenome_analysis",
        args=( 
            cell_id,
            jira_id,
            snv_jira_id,
            library_id,
            tantalus,
            prefix,
            out_dir,
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
    outprefix = os.path.join(out_dir,  cell_id + "_" + library_id)

    mutations_per_cell_png =  os.path.join(outprefix, "mutations_per_cell.png")
    summary_csv =  os.path.join(outprefix, "summary.csv")
    snvs_high_impact_csv =  os.path.join(outprefix, "snvs_high_impact.csv")
    snvs_all_csv =  os.path.join(outprefix, "snvs_all.csv")
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
            mgd.OutputFile(reporthtml),
            out_dir,
            cell_id + "_" + library_id,


        ),
    )

# def create_qc_workflow2(cell_id, jira_id, snv_jira_id, library_id, pseudobulk_groups, out_dir, reporthtml, mutationreporthtml, maf):
#     ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1, }

#     tantalus = "/work/shah/tantalus/"
#     prefixes = [os.path.join(out_dir, pseudobulk_group) for pseudobulk_group in pseudobulk_groups]

    # mutations_per_cell_pngs =  [os.path.join(prefix, "mutations_per_cell.png") for prefix in prefixes]

    # summary_csvs =  [os.path.join(prefix, "summary.csv") for prefix in prefixes]
    # snvs_high_impact_csvs =  [os.path.join(prefix, "snvs_high_impact.csv") for prefix in prefixes]
    # snvs_all_csvs =  [os.path.join(prefix, "snvs_all.csv") for prefix in prefixes]
    # trinuc_csvs =  [os.path.join(prefix, "trinuc.csv") for prefix in prefixes]
    # snv_adjacent_distance_pngs =  [os.path.join(prefix, "snv_adjacent_distance.png") for prefix in prefixes]
    # snv_genome_count_pngs =  [os.path.join(prefix,  "snv_genome_count.png")  for prefix in prefixes]
    # snv_cell_counts_pngs = [ os.path.join(prefix,  "snv_cell_counts.png") for prefix in prefixes]  
    # snv_alt_counts_pngs =  [os.path.join(prefix,  "snv_alt_counts.png")  for prefix in prefixes] 
    # rearranegementtype_distribution_pngs = [ os.path.join(prefix,  "rearranegementtype_distribution.png")   for prefix in prefixes]
    # chromosome_types_pngs =  [os.path.join(prefix,  "chromosome_types.png")  for prefix in prefixes]
    # BAFplot_pngs =  [os.path.join(prefix,  "BAFplot.png")  for prefix in prefixes]
    # CNplot_pngs = [ os.path.join(prefix,  "CNplot.png")  for prefix in prefixes]
    # datatype_summary_csvs =  [os.path.join(prefix,  "datatype_summary.csv") for prefix in prefixes]


    # workflow = pypeliner.workflow.Workflow(ctx=ctx)

    # workflow.transform(
    #     name='scgenome_plots',
    #     # ctx={'io': 1, 'mem': 8, 'disk': 100},
    #     func="single_cell.workflows.qc.scripts.scgenome-analysis.scgenome_analysis",
    #     args=( 

    #         cell_id,
    #         jira_id,
    #         snv_jira_id,
    #         library_id,
    #         tantalus,
    #         prefix,
    #         mgd.OutputFile(mutations_per_cell_png), 
    #         mgd.OutputFile(summary_csv), 
    #         mgd.OutputFile(snvs_high_impact_csv), 
    #         mgd.OutputFile(snvs_all_csv), 
    #         mgd.OutputFile(trinuc_csv), 
    #         mgd.OutputFile(snv_adjacent_distance_png), 
    #         mgd.OutputFile(snv_genome_count_png), 
    #         mgd.OutputFile(snv_cell_counts_png), 
    #         mgd.OutputFile(snv_alt_counts_png), 
    #         mgd.OutputFile(rearranegementtype_distribution_png), 
    #         mgd.OutputFile(chromosome_types_png),
    #         mgd.OutputFile(BAFplot_png), 
    #         mgd.OutputFile(CNplot_png), 
    #         mgd.OutputFile(datatype_summary_csv),

    #     ),
    # )

    # workflow.transform(
    #     name='create_main_report',
    #     # ctx={'io': 1, 'mem': 8, 'disk': 100},
    #     func="single_cell.workflows.qc.tasks.sample_level_report",
    #     args=(
    #         mgd.InputFile(mutations_per_cell_png), 
    #         mgd.InputFile(summary_csv), 
    #         mgd.InputFile(snvs_high_impact_csv), 
    #         mgd.InputFile(snvs_all_csv), 
    #         mgd.InputFile(trinuc_csv), 
    #         mgd.InputFile(snv_adjacent_distance_png), 
    #         mgd.InputFile(snv_genome_count_png), 
    #         mgd.InputFile(snv_cell_counts_png), 
    #         mgd.InputFile(snv_alt_counts_png), 
    #         mgd.InputFile(rearranegementtype_distribution_png), 
    #         mgd.InputFile(chromosome_types_png),
    #         mgd.InputFile(BAFplot_png), 
    #         mgd.InputFile(CNplot_png), 
    #         mgd.InputFile(datatype_summary_csv),


    #     ),
    # )
# rule createreport:
#     input:
#         snvaltcounts = "results/{pseudobulk_group}/{sample_id}_snv_alt_counts.png",
#         snvcellcounts = "results/{pseudobulk_group}/{sample_id}_snv_cell_counts.png",
#         trinuc = "results/{pseudobulk_group}/{sample_id}_trinuc.csv",
#         highimpact = "results/{pseudobulk_group}/{sample_id}_snvs_high_impact.csv",
#         allsnvs = "results/{pseudobulk_group}/{sample_id}_snvs_all.csv",
#         summary = "results/{pseudobulk_group}/{sample_id}_summary.csv",
#         datatype_summary = "results/{pseudobulk_group}/{sample_id}_datatype_summary.csv",
#         genome_count = "results/{pseudobulk_group}/{sample_id}_snv_genome_count.png",
#         adj_distance = "results/{pseudobulk_group}/{sample_id}_snv_adjacent_distance.png",
#         #breakpoint_adj_distance = "results/{pseudobulk_group}/{sample_id}_breakpoint_adjacent_distance.png",
#         breakpoint_dist = "results/{pseudobulk_group}/{sample_id}_rearranegementtype_distribution.png",
#         chr_types = "results/{pseudobulk_group}/{sample_id}_chromosome_types.png",
#         muts_per_cell = "results/{pseudobulk_group}/{sample_id}_mutations_per_cell.png",
#         BAFplot = "results/{pseudobulk_group}/{sample_id}_BAFplot.png",
#         CNplot = "results/{pseudobulk_group}/{sample_id}_CNplot.png",
#         indels = "results/indels/{pseudobulk_group}/{sample_id}.maf"
#     output: "reports/{pseudobulk_group}/{sample_id}.html"
#     params:
#         sample_id = "{sample_id}",
#         abs_snvaltcounts = lambda wildcards, input: os.path.abspath(input.snvaltcounts),
#         abs_snvcellcounts = lambda wildcards, input: os.path.abspath(input.snvcellcounts),
#         abs_trinuc = lambda wildcards, input: os.path.abspath(input.trinuc),
#         abs_highimpact = lambda wildcards, input: os.path.abspath(input.highimpact),
#         abs_allsnvs = lambda wildcards, input: os.path.abspath(input.allsnvs),
#         abs_summary = lambda wildcards, input: os.path.abspath(input.summary),
#         abs_datatype_summary = lambda wildcards, input: os.path.abspath(input.datatype_summary),
#         abs_genome_count = lambda wildcards, input: os.path.abspath(input.genome_count),
#         abs_muts_per_cell = lambda wildcards, input: os.path.abspath(input.muts_per_cell),
#         abs_adj_distance = lambda wildcards, input: os.path.abspath(input.adj_distance),
#         #abs_breapoint_adj_distance = lambda wildcards, input: os.path.abspath(input.breakpoint_adj_distance),
#         abs_chr_types = lambda wildcards, input: os.path.abspath(input.chr_types),
#         abs_breakpoint_dist = lambda wildcards, input: os.path.abspath(input.breakpoint_dist),
#         abs_BAFplot = lambda wildcards, input: os.path.abspath(input.BAFplot),
#         abs_CNplot = lambda wildcards, input: os.path.abspath(input.CNplot),
#         abs_indels = lambda wildcards, input: os.path.abspath(input.indels),
#     singularity:
#         "docker://rocker/tidyverse"
#     script:
#         "scripts/report.Rmd"


    return workflow


# def annotate_indel_vcfs(
#         workflow.transform(c
#         name='annotate_indels',
#         ctx=helpers.get_default_ctx(
#             memory=10,
#         ),
#         axesorigin
#         func="scripts.vcf2maf.py",
#         args=(
#             files=_get_indelvcfs()
#             vepdata,
#             genomeref,
#         kwargs={ "docker_image": vcf2maf}
#         )
#     )
# )

# def run_scgenome(
#         ticket, 
#         sample_id,
#         library_id,
#         snv_genotyping_ticket,
#         tantalus,
#         prefix,
#         trinuc,
#         snvs_high_impact,
#         snvs_all,
#         summary,
#         datatype_summary,
#         snv_genome_count,
#         snv_cell_counts,
#         snv_alt_counts,
#         BAFplot,
#         CNplot,
#         mutations_per_cell,
#         snv_adj_dist,
#         rearrangement_dist,
#         breakpoint_adj_dist,
#         chrom_types):


#     workflow.transform(
#         name='run_scgenome',
#         ctx=helpers.get_default_ctx(
#             memory=10,
#         ),
#         func="scripts.scgenome-analysis.py",
#         args=(
#             ticket, 
#             sample_id,
#             library_id,
#             snv_genotyping_ticket,
#             tantalus,
#             prefix,
#             mgd.OutputFile(trinuc),
#             mgd.OutputFile(snvs_high_impact),
#             mgd.OutputFile(snvs_all),
#             mgd.OutputFile(summary),
#             mgd.OutputFile(datatype_summary),
#             mgd.OutputFile(snv_genome_count),
#             mgd.OutputFile(snv_cell_counts),
#             mgd.OutputFile(snv_alt_counts),
#             mgd.OutputFile(BAFplot),
#             mgd.OutputFile(CNplot),
#             mgd.OutputFile(mutations_per_cell),
#             mgd.OutputFile(snv_adj_dist),
#             mgd.OutputFile(rearrangement_dist),
#             mgd.OutputFile(reakpoint_adj_dist),
#             mgd.OutputFile(chrom_types)
#         ),
#     )



# def postprocessing(        
#         sample_id,
#         single_node=False,
#         html_report,
#         output_maf,
#         mutation_report)

#     refdir_paths = config.refdir_data(refdir)['paths']
#     refdir_params = config.refdir_data(refdir)['params']

#     ideogram = refdir_paths["ideogram"]

#     titan_calls = titan[sample_id]
#     remixt_calls = remixt[sample_id]
#     sv_calls = breakpoints_consensus[sample_id]
#     roh_calls = roh[sample_id]
#     germline_vcf = germline_calls[sample_id]
#     somatic_calls = somatic_calls[sample_id]
#     chromosomes = refdir_params['chromosomes']

#     workflow = pypeliner.workflow.Workflow()

#     workflow.subworkflow(
#         name='coverage_normal_data',
#         func=get_coverage_data,
#         args=(
#             mgd.InputFile(normal_bam),x
#         ),
#         kwargs={'single_node': single_node}
#     )



#     workflow.subworkflow(
#         name='coverage_tumour_data',
#         func=get_coverage_data,
#         args=(
#             mgd.InputFile(tumour_bam),
#             mgd.TempOutputFile('tumour_coverage'),
#             refdir,
#         ),
#         kwargs={'single_node': single_node}
#     )


#     workflow.transform(
#         name='parse_roh',
#         ctx=helpers.get_default_ctx(
#             memory=5
#         ),
#         func="wgs.workflows.postprocessing.tasks.parse_roh",
#         args=(
#             mgd.InputFile(roh_calls),
#             mgd.TempOutputFile("ROH_parsed"),
#         ),
#     )

#     if remixt_calls:

#         workflow.transform(
#             name='generate_genome_wide_plot',
#             ctx=helpers.get_default_ctx(
#                 memory=10,
#             ),
#             func="wgs.workflows.postprocessing.tasks.genome_wide",
#             args=(
#                 mgd.InputFile(titan_calls),
#                 mgd.TempInputFile("ROH_parsed"),
#                 mgd.InputFile(germline_vcf),
#                 mgd.InputFile(somatic_calls),
#                 mgd.TempInputFile('tumour_coverage'),
#                 mgd.TempInputFile('normal_coverage'),
#                 mgd.InputFile(sv_calls),
#                 mgd.InputFile(ideogram),
#                 chromosomes,
#                 mgd.OutputFile(genome_wide_plot),

#             ),
#             kwargs={"remixt": mgd.InputFile(remixt_calls),
#                     "remixt_label": sample_id}
#         )
#         workflow.transform(
#             name='generate_circos_plot',
#             ctx=helpers.get_default_ctx(
#                 memory=10
#             ),
#             func="wgs.workflows.postprocessing.tasks.circos",
#             args=(
#                 mgd.InputFile(titan_calls),
#                 sample_id,
#                 mgd.InputFile(sv_calls),
#                 mgd.TempOutputFile(circos_plot_remixt),
#                 mgd.TempOutputFile(circos_plot_titan),
#                 mgd.TempSpace('circos'),
#             ),

#             kwargs={'docker_image': config.containers('circos'),
#                     'remixt_calls': mgd.InputFile(remixt_calls)},
#         )
#     else:

#         workflow.transform(
#             name='generate_genome_wide_plot',
#             ctx=helpers.get_default_ctx(
#                 memory=10,
#             ),
#             func="wgs.workflows.postprocessing.tasks.genome_wide",
#             args=(
#                 mgd.InputFile(titan_calls),
#                 mgd.TempInputFile("ROH_parsed"),
#                 mgd.InputFile(germline_vcf),
#                 mgd.InputFile(somatic_calls),
#                 mgd.TempInputFile('tumour_coverage'),
#                 mgd.TempInputFile('normal_coverage'),
#                 mgd.InputFile(sv_calls),
#                 mgd.InputFile(ideogram),
#                 chromosomes,
#                 mgd.OutputFile(genome_wide_plot),
#             ),
#         )

#         workflow.transform(
#             name='generate_circos_plot',
#             ctx=helpers.get_default_ctx(
#                 memory=10
#             ),
#             func="wgs.workflows.postprocessing.tasks.circos",
#             args=(
#                 mgd.InputFile(titan_calls),
#                 sample_id,
#                 mgd.InputFile(sv_calls),
#                 mgd.TempOutputFile(circos_plot_remixt),
#                 mgd.TempOutputFile(circos_plot_titan),
#                 mgd.TempSpace('circos'),
#             ),

#             kwargs={'docker_image': config.containers('circos')}
#         )

#     return workflow



# df = pd.read_csv(config["metadata"])
# df['id'] = df['sample_id'] + '_' + df['library_id']
# df.set_index('id', inplace=True, drop=False)

# def _get_indelvcfs(wildcards):
#     jiraid = df.loc[wildcards.sample_id, 'jira_id']
#     dir1 = config["tantalus"] + jiraid + "/results/variants/" + "*_strelka_indel.vcf.gz"
#     files1 = sorted(glob.glob(dir1))
#     dir2 = config["tantalus"] + jiraid + "/results/" + "*_strelka_indel.vcf.gz"
#     files2 = sorted(glob.glob(dir2))
#     dir3 = config["tantalus"] + jiraid + "/results/variant_calling/*/" + "*indel.vcf.gz"
#     files3 = sorted(glob.glob(dir3))
#     files = files1 + files2 + files3
#     return files

# rule all:
#     input:
#         expand("reports/{pseudobulk_group}/{sample_id}.html", zip, pseudobulk_group=df["pseudobulk_group"], sample_id=df['id']),
#         expand("results/indels/{pseudobulk_group}/{sample_id}.maf", zip, pseudobulk_group=df["pseudobulk_group"], sample_id=df['id']),
#         expand("mutationreports/{pseudobulk_group}.html", pseudobulk_group=df["pseudobulk_group"].unique())



# # rule runscgenome:
# #     output:
# #         "results/{pseudobulk_group}/{sample_id}_trinuc.csv",
# #         "results/{pseudobulk_group}/{sample_id}_snvs_high_impact.csv",
# #         "results/{pseudobulk_group}/{sample_id}_snvs_all.csv",
# #         "results/{pseudobulk_group}/{sample_id}_summary.csv",
# #         "results/{pseudobulk_group}/{sample_id}_datatype_summary.csv",
# #         "results/{pseudobulk_group}/{sample_id}_snv_genome_count.png",
# #         "results/{pseudobulk_group}/{sample_id}_snv_cell_counts.png",
# #         "results/{pseudobulk_group}/{sample_id}_snv_alt_counts.png",
# #         "results/{pseudobulk_group}/{sample_id}_BAFplot.png",
# #         "results/{pseudobulk_group}/{sample_id}_CNplot.png",
# #         "results/{pseudobulk_group}/{sample_id}_mutations_per_cell.png",
# #         "results/{pseudobulk_group}/{sample_id}_snv_adjacent_distance.png",
# #         "results/{pseudobulk_group}/{sample_id}_rearranegementtype_distribution.png",
# #         #"results/{pseudobulk_group}/{sample_id}_breakpoint_adjacent_distance.png",
# #         "results/{pseudobulk_group}/{sample_id}_chromosome_types.png",
# #     params:
# #         ticket = lambda wildcards: df.loc[wildcards.sample_id, 'jira_id'],
# #         sample_id = lambda wildcards: df.loc[wildcards.sample_id, 'sample_id'],
# #         library_id = lambda wildcards: df.loc[wildcards.sample_id, 'library_id'],
# #         snv_genotyping_ticket = lambda wildcards: df.loc[wildcards.sample_id, 'snv_genotyping_jira_id'],
# #         tantalus = config["tantalus"],
# #         prefix = "results/{pseudobulk_group}/{sample_id}"
# #     script:
# #         "scripts/scgenome-analysis.py"


# # checkpoint annotateindels:
# #     input: _get_indelvcfs
# #     output: directory("results/indelstemp/{pseudobulk_group}/{sample_id}/maf/")
# #     params:
# #         vepdata = config["vepdata"],
# #         genomeref = config["genomeref"]
# #     singularity:
# #         "shub://rdmorin/cancer_docker_singularity:vcf2maf"
# #     script: "scripts/vcf2maf.py"

# def aggregate_input(wildcards):
#     checkpoint_output = checkpoints.annotateindels.get(**wildcards).output[0]
#     maffiles = expand("results/indelstemp/{pseudobulk_group}/{sample_id}/maf/{sample_id}_{indelid}.maf",
#            pseudobulk_group=wildcards.pseudobulk_group,
#            sample_id = wildcards.sample_id,
#            indelid=glob_wildcards(os.path.join(checkpoint_output, "{samp}_{indelid}.maf")).indelid)
#     return maffiles

# def run_scgenome(
#         ticket, 
#         sample_id,
#         library_id,
#         snv_genotyping_ticket,
#         tantalus,
#         prefix,
#         trinuc,
#         snvs_high_impact,
#         snvs_all,
#         summary,
#         datatype_summary,
#         snv_genome_count,
#         snv_cell_counts,
#         snv_alt_counts,
#         BAFplot,
#         CNplot,
#         mutations_per_cell,
#         snv_adj_dist,
#         rearrangement_dist,
#         breakpoint_adj_dist,
#         chrom_types):


#     workflow.transform(
#         name='merge_indels',
#         ctx=helpers.get_default_ctx(
#             memory=10,
#         ),
#         axesoriginrev
#         func="scripts.scgenome-analysis.py",
#         args=(
#             files=aggregate_input,
#         ),
#     )


# rule mergeindels:
#     input: aggregate_input
#     output: temp("results/indels/{pseudobulk_group}/{sample_id}.maf")
#     shell:
#         """
#         awk '(NR == 2) || (FNR > 2)' {input} > {output}
#         """

# def run_scgenome(
#         ticket, 
#         sample_id,
#         library_id,
#         snv_genotyping_ticket,
#         tantalus,
#         prefix,
#         trinuc,
#         snvs_high_impact,
#         snvs_all,
#         summary,
#         datatype_summary,
#         snv_genome_count,
#         snv_cell_counts,
#         snv_alt_counts,
#         BAFplot,
#         CNplot,
#         mutations_per_cell,
#         snv_adj_dist,
#         rearrangement_dist,
#         breakpoint_adj_dist,
#         chrom_types):


#     workflow.transform(
#         name='create_report',
#         ctx=helpers.get_default_ctx(
#             memory=10,
#         ),
#         func="scripts.report.Rmd",
#         args=(
#             sample_id,
#             abs_snvaltcounts,
#             abs_snvcellcounts,
#             abs_trinuc ,
#             abs_highimpact, 
#             abs_allsnvs ,
#             abs_summary ,
#             abs_datatype_summary ,
#             abs_genome_count ,
#             abs_muts_per_cell,
#             abs_adj_distance ,
#             abs_chr_types,
#             abs_breakpoint_dist, 
#             abs_BAFplot,
#             abs_CNplot,
#             abs_indels 
#         ),
        
#         kwargs={"docker_image": tidyverse}
#     )

# rule createreport:
#     input:
#         snvaltcounts = "results/{pseudobulk_group}/{sample_id}_snv_alt_counts.png", 
#         snvcellcounts = "results/{pseudobulk_group}/{sample_id}_snv_cell_counts.png",
#         trinuc = "results/{pseudobulk_group}/{sample_id}_trinuc.csv",
#         highimpact = "results/{pseudobulk_group}/{sample_id}_snvs_high_impact.csv",
#         allsnvs = "results/{pseudobulk_group}/{sample_id}_snvs_all.csv",
#         summary = "results/{pseudobulk_group}/{sample_id}_summary.csv",
#         datatype_summary = "results/{pseudobulk_group}/{sample_id}_datatype_summary.csv",
#         genome_count = "results/{pseudobulk_group}/{sample_id}_snv_genome_count.png",
#         adj_distance = "results/{pseudobulk_group}/{sample_id}_snv_adjacent_distance.png",
#         #breakpoint_adj_distance = "results/{pseudobulk_group}/{sample_id}_breakpoint_adjacent_distance.png",
#         breakpoint_dist = "results/{pseudobulk_group}/{sample_id}_rearranegementtype_distribution.png",
#         chr_types = "results/{pseudobulk_group}/{sample_id}_chromosome_types.png",
#         muts_per_cell = "results/{pseudobulk_group}/{sample_id}_mutations_per_cell.png",
#         BAFplot = "results/{pseudobulk_group}/{sample_id}_BAFplot.png",
#         CNplot = "results/{pseudobulk_group}/{sample_id}_CNplot.png",
#         indels = "results/indels/{pseudobulk_group}/{sample_id}.maf"
#     output: "reports/{pseudobulk_group}/{sample_id}.html"
#     params:
#         sample_id = "{sample_id}",
#         abs_snvaltcounts = lambda wildcards, input: os.path.abspath(input.snvaltcounts),
#         abs_snvcellcounts = lambda wildcards, input: os.path.abspath(input.snvcellcounts),
#         abs_trinuc = lambda wildcards, input: os.path.abspath(input.trinuc),
#         abs_highimpact = lambda wildcards, input: os.path.abspath(input.highimpact),
#         abs_allsnvs = lambda wildcards, input: os.path.abspath(input.allsnvs),
#         abs_summary = lambda wildcards, input: os.path.abspath(input.summary),
#         abs_datatype_summary = lambda wildcards, input: os.path.abspath(input.datatype_summary),
#         abs_genome_count = lambda wildcards, input: os.path.abspath(input.genome_count),
#         abs_muts_per_cell = lambda wildcards, input: os.path.abspath(input.muts_per_cell),
#         abs_adj_distance = lambda wildcards, input: os.path.abspath(input.adj_distance),
#         #abs_breapoint_adj_distance = lambda wildcards, input: os.path.abspath(input.breakpoint_adj_distance),
#         abs_chr_types = lambda wildcards, input: os.path.abspath(input.chr_types),
#         abs_breakpoint_dist = lambda wildcards, input: os.path.abspath(input.breakpoint_dist),
#         abs_BAFplot = lambda wildcards, input: os.path.abspath(input.BAFplot),
#         abs_CNplot = lambda wildcards, input: os.path.abspath(input.CNplot),
#         abs_indels = lambda wildcards, input: os.path.abspath(input.indels),
#     singularity:
#         "docker://rocker/tidyverse"
#     script:
#         "scripts/report.Rmd"


# def getindelfiles(wildcards):
#     files = expand("results/indels/{pseudobulk_group}/{sample}.maf", pseudobulk_group = wildcards.pseudobulk_group, sample = df[df["pseudobulk_group"] == wildcards.pseudobulk_group].index)
#     return files

# rule perpseudobulk_indelfile:
#     input: getindelfiles
#     output:
#         "results/indels_merged/maffiles/{pseudobulk_group}.maf",
#         "results/indels_merged/maffiles/{pseudobulk_group}-highimpact.maf"
#     singularity:
#         "docker://rocker/tidyverse"
#     script:
#         "scripts/mergemafs.R"

# rule createindelreport:
#     input:
#         indels="results/indels_merged/maffiles/{pseudobulk_group}.maf",
#         indels_highimpact="results/indels_merged/maffiles/{pseudobulk_group}-highimpact.maf"
#     params:
#         abs_indels = lambda wildcards, input: os.path.abspath(input.indels),
#         abs_indels_highimpact = lambda wildcards, input: os.path.abspath(input.indels_highimpact),
#     output:
#         "indelreports/{pseudobulk_group}.html"
#     singularity:
#         "docker://rocker/tidyverse"
#     script:
#         "scripts/indelreport.Rmd"

# def getsnvfiles(wildcards):
#     files = expand("results/{pseudobulk_group}/{sample}_snvs_all.csv", pseudobulk_group = wildcards.pseudobulk_group, sample = df[df["pseudobulk_group"] == wildcards.pseudobulk_group].index)
#     return files

# rule perpseudobulk_snvfile:
#     input: getsnvfiles
#     output:
#         "results/snvs_merged/csvfiles/{pseudobulk_group}.csv",
#     singularity:
#         "docker://rocker/tidyverse"
#     script:
#         "scripts/mergesnvs.R"

# rule createmutationreport:
#     input:
#         snvs="results/snvs_merged/csvfiles/{pseudobulk_group}.csv",
#         indels="results/indels_merged/maffiles/{pseudobulk_group}.maf",
#         indels_highimpact="results/indels_merged/maffiles/{pseudobulk_group}-highimpact.maf"
#     params:
#         abs_snvs = lambda wildcards, input: os.path.abspath(input.snvs),
#         abs_indels = lambda wildcards, input: os.path.abspath(input.indels),
#         abs_indels_highimpact = lambda wildcards, input: os.path.abspath(input.indels_highimpact),
#     output:
#         "mutationreports/{pseudobulk_group}.html"
#     singularity:
#         "docker://rocker/tidyverse"
#     script:
#         "scripts/mutationreport.Rmd"

# annotateindels
# mergeindels
# createreport
# perpseudobulk_indelfile
# createindelreport
# perpseudobulk_snvfile
# createmutationreport

