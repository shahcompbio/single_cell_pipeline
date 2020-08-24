import os

import pandas as pd
from single_cell.utils import helpers

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


def filter_snvs_for_high_impact(snv, filtsnv, docker_image=None):
    cmd = ["mergesnvs.R", snv, filtsnv]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def filter_maf_for_high_impact(maf, filtmaf, docker_image=None):
    cmd = ["mergemafs.R", maf, filtmaf]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def vcf2maf(vcf_file, output_maf, tempdir, vep_ref, docker_image=None):
    if vcf_file.endswith('.gz'):
        helpers.makedirs(tempdir)
        vcf_unzipped = os.path.join(tempdir, 'unzipped_vcf.vcf')
        helpers.gunzip_file(vcf_file, vcf_unzipped)
    else:
        vcf_unzipped = vcf_file

    cmd = [
        'vcf2maf.pl', '--input-vcf', vcf_unzipped, '--output-maf', output_maf,
        '--vep-path', '/usr/local/bin',
        '--ref-fasta', vep_ref['reference_fasta'],
        '--filter-vcf', vep_ref['reference_filter_vcf'],
        '--vep-data', vep_ref['reference_dir'],
    ]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def sample_level_report(
        mutations_per_cell, summary,
        snvs_high_impact, snvs_all, trinuc, snv_adjacent_distance, snv_genome_count,
        snv_cell_counts, snv_alt_counts, rearranegementtype_distribution_destruct_unfiltered,
        chromosome_types_destruct_unfiltered, rearranegementtype_distribution_destruct_filtered,
        chromosome_types_destruct_filtered, rearranegementtype_distribution_lumpy_unfiltered,
        chromosome_types_lumpy_unfiltered, baf_plot, cn_plot, datatype_summary, maf, html_file,
        sample_id, docker_image=None
):
    cmd = [
        'run_report.sh', html_file, sample_id, mutations_per_cell, summary,
        snvs_high_impact, snvs_all, trinuc, snv_adjacent_distance, snv_genome_count,
        snv_cell_counts, snv_alt_counts, rearranegementtype_distribution_destruct_unfiltered,
        chromosome_types_destruct_unfiltered, rearranegementtype_distribution_destruct_filtered,
        chromosome_types_destruct_filtered, rearranegementtype_distribution_lumpy_unfiltered,
        chromosome_types_lumpy_unfiltered, baf_plot, cn_plot, datatype_summary, maf
    ]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def create_mutation_report(
        pseudo_bulk_group, merged_maf, high_impact_maf, high_impact_snvs, report_html,
        docker_image=None
):
    cmd = [
        "run_mutationreport.sh", report_html, pseudo_bulk_group,
        high_impact_snvs, merged_maf, high_impact_maf
    ]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)
