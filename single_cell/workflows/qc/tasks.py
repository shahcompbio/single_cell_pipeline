import os

import pandas as pd
from single_cell.utils import helpers
import tarfile
import pypeliner.commandline


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


def filter_snvs_for_high_impact(snv, filtsnv):
    cmd = ["mergesnvs.R", snv, filtsnv]
    pypeliner.commandline.execute(*cmd)


def filter_maf_for_high_impact(maf, filtmaf):
    cmd = ["mergemafs.R", maf, filtmaf]
    pypeliner.commandline.execute(*cmd)


def vcf2maf(vcf_file, output_maf, tempdir, vep_ref, docker_image):
    if vcf_file.endswith('.gz'):
        helpers.makedirs(tempdir)
        vcf_unzipped = os.path.join(tempdir, 'unzipped_vcf.vcf')
        helpers.gunzip_file(vcf_file, vcf_unzipped)
    else:
        vcf_unzipped = vcf_file

    # cmd = [
    #     'vcf2maf.pl', '--input-vcf', vcf_unzipped, '--output-maf', output_maf,
    #     '--vep-path', '/usr/local/bin',
    #     '--ref-fasta', vep_ref['reference_fasta'],
    #     '--filter-vcf', vep_ref['reference_filter_vcf'],
    #     '--vep-data', vep_ref['reference_dir'],
    vepdata = "/work/shah/reference/vep/"
    # cmd = [  'vcf2maf.pl', '--input-vcf', vcf_unzipped, '--output-maf', output_maf,
    #     '--vep-path', '/home/abramsd/miniconda3/envs/r-environment/bin',
    #     '--ref-fasta',
    #     os.path.join(reference, 'homo_sapiens', '99_GRCh37', 'Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz'),
    #     '--filter-vcf', os.path.join(reference, 'ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz'),
    #     '--vep-data', reference,
    # ]
    reference = "/juno/work/shah/reference/genomes/GRCh37-lite/GRCh37-lite.fa"
    cmd = ["vcf2maf.pl", "--input-vcf", vcf_unzipped, "--output-maf", output_maf,
           "--ref-fasta", reference, "--filter-vcf", "0", "--vep-path", "/home/abramsd/miniconda3/envs/r-environment/bin", 
           "--vep-data", vepdata]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def sample_level_report(
    snvs_all, plots_tar, maf, html_file, sample_id
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

    cmd = [
        'run_report.sh', *list(map(os.path.abspath, [html_file, sample_id, mutations_per_cell, summary,
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
        baf_plot, cn_plot, datatype_summary, maf]))
    ]
    
    pypeliner.commandline.execute(*cmd)
    tar.close()    


def create_mutation_report(
        pseudobulk_group, merged_maf, high_impact_maf, high_impact_snvs, report_html
):
    cmd = [
        "run_mutationreport.sh", *list(map(os.path.abspath, [report_html, pseudobulk_group,
        high_impact_snvs, merged_maf, high_impact_maf]))
    ]
    pypeliner.commandline.execute(*cmd)
