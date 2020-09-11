import os
import shutil

import pandas as pd
import pypeliner.commandline
from single_cell.utils import helpers
import tarfile


def merge_mafs(mafs, merged_maf, id_colname=False):
    assert isinstance(mafs, dict)
    preppedmafs = []
    for samplegroup, maf in mafs.items():
        maftemp = pd.read_csv(maf, sep="\t", skiprows=1)
        maftemp = maftemp.astype(
            {"Chromosome": "str", "PHENO": "str", "SOMATIC": "str", "PUBMED": "str"}
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


def filter_snvs_for_high_impact(snv, filtsnv, docker_image=None):
    cmd = ["mergesnvs.R", snv, filtsnv]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def filter_maf_for_high_impact(maf, filtmaf, docker_image=None):
    cmd = ["mergemafs.R", maf, filtmaf]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def vcf2maf(vcf_file, output_maf, tempdir, vep_ref, docker_image=None):
    vcf_file_copy = os.path.join(tempdir, 'vcf2maf_input_vcf', os.path.basename(vcf_file))
    helpers.makedirs(vcf_file_copy, isfile=True)
    shutil.copyfile(vcf_file, vcf_file_copy)

    if vcf_file_copy.endswith('.gz'):
        vcf_unzipped = os.path.join(tempdir, 'unzipped_vcf.vcf')
        helpers.gunzip_file(vcf_file_copy, vcf_unzipped)
    else:
        vcf_unzipped = vcf_file_copy

    # cmd = [
    #     'vcf2maf', vcf_unzipped, output_maf,
    #     vep_ref['reference_fasta'],
    #     vep_ref['reference_filter_vcf'],
    #     vep_ref['reference_dir'],
    # ]
    cmd = ["vcf2maf.pl", "--input-vcf", vcf_unzipped, "--output-maf", output_maf,
           "--ref-fasta", "/juno/work/shah/reference/genomes/GRCh37-lite/GRCh37-lite.fa", "--filter-vcf", "0", "--vep-path", "/home/abramsd/miniconda3/envs/r-environment/bin",
           "--vep-data", "/work/shah/reference/vep/"]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def sample_level_report(
    snvs_all, plots_tar, maf, html_file, sample_id, docker_image=None
):
    
    basedir = os.path.dirname(plots_tar)

    tar = tarfile.open(plots_tar)
    tar.extractall(basedir)
    

    type_distribution_destruct_unfiltspace = os.path.join(basedir, "qc_plots", 
        "destruct_unfilt", "type_distribution.png")
    rearrangement_type_distribution_destruct_unfiltspace = os.path.join(basedir, "qc_plots", 
        "destruct_unfilt",  "rearrangement_type_distribution.png")
    type_genome_distribution_destruct_unfiltspace = os.path.join(basedir, "qc_plots", 
        "destruct_unfilt",  "type_genome_distribution.png")
    rearrangement_type_genome_distribution_destruct_unfiltspace = os.path.join(basedir, "qc_plots", 
        "destruct_unfilt",  "rearrangement_type_genome_distribution.png")
    type_size_distribution_destruct_unfiltspace = os.path.join(basedir, "qc_plots", 
        "destruct_unfilt",  "type_size_distribution.png")
    rearrangement_type_size_distribution_destruct_unfiltspace = os.path.join(basedir, "qc_plots", 
        "destruct_unfilt",  "rearrangement_type_size_distribution.png")
    type_distribution_destruct_filtspace = os.path.join(basedir, "qc_plots", 
        "destruct_filt",  "type_distribution.png")
    rearrangement_type_distribution_destruct_filtspace = os.path.join(basedir, "qc_plots", 
        "destruct_filt",  "rearrangement_type_distribution.png")
    type_genome_distribution_destruct_filtspace = os.path.join(basedir, "qc_plots", 
        "destruct_filt",  "type_genome_distribution.png")
    rearrangement_type_genome_distribution_destruct_filtspace = os.path.join(basedir, "qc_plots", 
        "destruct_filt",  "rearrangement_type_genome_distribution.png")
    type_size_distribution_destruct_filtspace = os.path.join(basedir, "qc_plots", 
        "destruct_filt",  "type_size_distribution.png")
    rearrangement_type_size_distribution_destruct_filtspace = os.path.join(basedir, "qc_plots", 
        "destruct_filt",  "rearrangement_type_size_distribution.png")

    type_distribution_lumpy_unfiltspace = os.path.join(basedir, "qc_plots", 
        "lumpy_unfilt",  "type_distribution.png")
    rearrangement_type_distribution_lumpy_unfiltspace = os.path.join(basedir, "qc_plots", 
        "lumpy_unfilt",  "rearrangement_type_distribution.png")
    type_genome_distribution_lumpy_unfiltspace = os.path.join(basedir, "qc_plots", 
        "lumpy_unfilt", "type_genome_distribution.png")
    rearrangement_type_genome_distribution_lumpy_unfiltspace = os.path.join(basedir, 
        "qc_plots", "lumpy_unfilt", "rearrangement_type_genome_distribution.png")
    type_size_distribution_lumpy_unfiltspace = os.path.join(basedir, "qc_plots", 
        "lumpy_unfilt",  "type_size_distribution.png")
    rearrangement_type_size_distribution_lumpy_unfiltspace = os.path.join(basedir, "qc_plots", 
        "lumpy_unfilt",  "rearrangement_type_size_distribution.png")

    
    mutations_per_cell = os.path.join(basedir, "qc_plots", "mutations_per_cell.png")
    summary = os.path.join(basedir, "qc_plots", "summary.csv")
    snvs_high_impact = os.path.join(basedir, "qc_plots", "snvs_high_impact.csv")
    trinuc = os.path.join(basedir, "qc_plots", "trinuc.csv")
    snv_adjacent_distance = os.path.join(basedir, "qc_plots", "snv_adjacent_distance.png")
    snv_genome_count = os.path.join(basedir, "qc_plots", "snv_genome_count.png")
    snv_cell_counts = os.path.join(basedir, "qc_plots", "snv_cell_counts.png")
    snv_alt_counts = os.path.join(basedir, "qc_plots", "snv_alt_counts.png")
    baf_plot = os.path.join(basedir, "qc_plots", "baf_plot.png")
    cn_plot = os.path.join(basedir, "qc_plots", "cn_plot.png")
    datatype_summary = os.path.join(basedir, "qc_plots", "datatype_summary.csv")

    files_args = [html_file, sample_id, mutations_per_cell, summary,
        snvs_high_impact, snvs_all, trinuc, snv_adjacent_distance, snv_genome_count,
        snv_cell_counts, snv_alt_counts, type_distribution_destruct_unfiltspace, 
        rearrangement_type_distribution_destruct_unfiltspace, type_genome_distribution_destruct_unfiltspace,
        rearrangement_type_genome_distribution_destruct_unfiltspace, type_size_distribution_destruct_unfiltspace,
        rearrangement_type_size_distribution_destruct_unfiltspace, type_distribution_destruct_filtspace, rearrangement_type_distribution_destruct_filtspace,
        type_genome_distribution_destruct_filtspace,  rearrangement_type_genome_distribution_destruct_filtspace,type_size_distribution_destruct_filtspace, 
        rearrangement_type_size_distribution_destruct_filtspace, type_distribution_lumpy_unfiltspace, 
        rearrangement_type_distribution_lumpy_unfiltspace, type_genome_distribution_lumpy_unfiltspace, 
        rearrangement_type_genome_distribution_lumpy_unfiltspace, 
        type_size_distribution_lumpy_unfiltspace, rearrangement_type_size_distribution_lumpy_unfiltspace, 
        baf_plot, cn_plot, datatype_summary, maf
    ]

    files_args = [os.path.abspath(v) for v in files_args]

    cmd = ['run_report.sh'] + files_args 

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    tar.close()  


def create_mutation_report(
        pseudo_bulk_group, merged_maf, high_impact_maf, high_impact_snvs, report_html,
        docker_image=None
):
    cmd = [
        "run_mutationreport.sh",
        os.path.abspath(report_html),
        os.path.abspath(pseudo_bulk_group),
        os.path.abspath(high_impact_snvs),
        os.path.abspath(merged_maf),
        os.path.abspath(high_impact_maf)
    ]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)
