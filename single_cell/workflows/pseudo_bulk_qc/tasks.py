import os
import shutil
import pandas as pd
from single_cell.utils import helpers
import pypeliner.commandline
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
rmarkdown = importr("rmarkdown")

def merge_mafs(mafs, merged_maf, id_colname=False):
    assert isinstance(mafs, dict)
    preppedmafs = []

    for samplegroup, maf in mafs.items():
        maftemp = pd.read_csv(maf, sep="\t", skiprows=1, dtype='str')
        maftemp = maftemp.astype(
            {"Chromosome": "str", "PHENO": "str", "SOMATIC": "str", "PUBMED": "str", "t_alt_count": "int64"}
        )
        if id_colname:
            maftemp["id"] = [samplegroup] * len(maftemp)
        preppedmafs.append(maftemp)

    output = pd.concat(preppedmafs).reset_index()
    output = output[output.t_alt_count > 2]
    output.to_csv(merged_maf, sep="\t", index=False, header=True)


def merge_snvs(snv_files, merged_snv, id_colname=False):
    assert isinstance(snv_files, dict)

    preppednsvs = []
    for samplegroup, snv_file in snv_files.items():
        svtemp = pd.read_csv(snv_file).astype({"chrom": "str"})
        if id_colname:
            svtemp["id"] = [samplegroup] * len(svtemp)
        svtemp = svtemp[svtemp.alt_counts > 1]
        preppednsvs.append(svtemp)

    output = pd.concat(preppednsvs).reset_index()
    output.to_csv(merged_snv, sep="\t", index=False, header=True)


def filter_snvs_for_high_impact(snv, filtsnv):
    script_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'scripts','mergesnvs.R')
    cmd = ["Rscript", script_path, snv, filtsnv]
    pypeliner.commandline.execute(*cmd)


def filter_maf_for_high_impact(maf, filtmaf):
    script_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'scripts','mergemafs.R')
    cmd = ["Rscript", script_path, maf, filtmaf]
    pypeliner.commandline.execute(*cmd)


def vcf2maf(vcf_file, output_maf, tempdir, vep_ref):
    vcf_file_copy = os.path.join(tempdir, 'vcf2maf_input_vcf', os.path.basename(vcf_file))
    helpers.makedirs(vcf_file_copy, isfile=True)
    shutil.copyfile(vcf_file, vcf_file_copy)

    if vcf_file_copy.endswith('.gz'):
        vcf_unzipped = os.path.join(tempdir, 'unzipped_vcf.vcf')
        helpers.gunzip_file(vcf_file_copy, vcf_unzipped)
    else:
        vcf_unzipped = vcf_file_copy

    cmd = [
        'vcf2maf', 
        vcf_unzipped, 
        output_maf,
        vep_ref['reference_fasta'],
        vep_ref['reference_filter_vcf'],
        vep_ref['reference_dir'],
    ]
    
    pypeliner.commandline.execute(*cmd)


def sample_level_report(
        mutations_per_cell, summary, snvs_high_impact, snvs_all,
        trinuc, snv_adjacent_distance, snv_genome_count,
        snv_cell_counts, snv_alt_counts, destruct_rearrangement_plots_unfiltered,
        destruct_rearrangement_plots_filtered, lumpy_rearrangement_plots_unfiltered,
        baf_plot, cn_plot, datatype_summary, maf, html_file,
        sample_id, tmpspace
):

    rmd_script = os.path.join(os.path.dirname(os.path.realpath(__file__)), 
        'scripts','report.Rmd'
    )

    parameters = robjects.r.list(sample_id=sample_id, mutations_per_cell_png=mutations_per_cell,
        summary_csv=summary, snvs_high_impact_csv=snvs_high_impact, snvs_all_csv=snvs_all, 
        trinuc_csv=trinuc,  snv_adjacent_distance_png=snv_adjacent_distance, snv_genome_count_png=snv_genome_count, 
        snv_cell_counts_png=snv_cell_counts, snv_alt_counts_png=snv_alt_counts, 
        destruct_rearrangement_plots_unfiltered=destruct_rearrangement_plots_unfiltered,
        destruct_rearrangement_plots_filtered=destruct_rearrangement_plots_filtered, 
        lumpy_rearrangement_plots_unfiltered=lumpy_rearrangement_plots_unfiltered,
        BAFplot_png=baf_plot, cn_plot_png=cn_plot, datatype_summary_csv=datatype_summary, 
        maf=maf
    )

    rmarkdown.render(rmd_script, output_file=html_file, 
        output_options=robjects.r.list(self_contained=True), params=parameters
    )
    

def create_mutation_report(
        pseudo_bulk_group, merged_maf, high_impact_maf,
        high_impact_snvs, report_html, tmpspace
):
    rmd_script = os.path.join(os.path.dirname(os.path.realpath(__file__)), 
        'scripts','mutationreport.Rmd'
    )

    parameters = robjects.r.list(pseudobulk_group=pseudo_bulk_group,
        merged_filt_snvs=high_impact_snvs, merged_maf=merged_maf, 
        high_impact_maf=high_impact_maf
    )

    rmarkdown.render(rmd_script, output_file=report_html, 
        output_options=robjects.r.list(self_contained=True), params=parameters
    )
    