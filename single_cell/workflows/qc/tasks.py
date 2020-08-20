import os
import pypeliner.commandline 
import glob
import gzip
import shutil
from single_cell.utils import helpers
import pandas as pd


def merge_mafs(mafs, merged_maf, id_colname=False):
    assert isinstance(mafs, dict)
    preppedmafs = []
    for samplegroup, maf in mafs.items():
        maftemp = pd.read_csv(maf, sep="\t", skiprows=1)
        maftemp = maftemp.astype({"Chromosome": "str","PHENO": "str","SOMATIC": "str","PUBMED": "str"})
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
            svtemp["id"] = [samplegroup] *len(svtemp)
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


def _get_indelvcfs(dir):
    indel_vcfs = []
    indel_vcfs += sorted(glob.glob(dir + "/results/variants/" + "*_strelka_indel.vcf.gz" ))
    indel_vcfs += sorted(glob.glob(dir + "/results/" + "*_strelka_indel.vcf.gz" ))
    indel_vcfs += sorted(glob.glob(dir + "/results/variant_calling/*/" + "*indel.vcf.gz" ))
    return indel_vcfs


def get_snv_all_csvs(dir):
    output = {}
    dir_len = len(dir.split("/"))
    filenames =  sorted(glob.glob(dir + "/*/*/*snv*all*csv", recursive=True))
    for filename in filenames:

        name_len = len(filename.split("/"))
        samplelabel = filename.split("/")[dir_len : name_len - 1]
        output[tuple(samplelabel)] = filename
    return output


def vcf2maf(vcf_file, output_maf, tempdir, reference, docker_image):

    if vcf_file.endswith('.gz'):
        helpers.makedirs(tempdir)
        vcf_unzipped = os.path.join(tempdir, 'unzipped_vcf.vcf')
        helpers.gunzip_file(vcf_file, vcf_unzipped)
    else:
        vcf_unzipped = vcf_file

    cmd = [  'vcf2maf.pl', '--input-vcf', vcf_unzipped, '--output-maf', output_maf,
        '--vep-path', '/usr/local/bin',
        '--ref-fasta',
        os.path.join(reference, 'homo_sapiens', '99_GRCh37', 'Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz'),
        '--filter-vcf', os.path.join(reference, 'ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz'),
        '--vep-data', reference,
    ]

    pypeliner.commandline.execute(*cmd, docker_image = docker_image)


def sample_level_report(mutations_per_cell, summary, 
    snvs_high_impact, snvs_all, trinuc, snv_adjacent_distance, snv_genome_count, 
    snv_cell_counts, snv_alt_counts, rearranegementtype_distribution_destruct_unfiltered, 
    chromosome_types_destruct_unfiltered,rearranegementtype_distribution_destruct_filtered, 
    chromosome_types_destruct_filtered,rearranegementtype_distribution_lumpy_unfiltered, 
    chromosome_types_lumpy_unfiltered,
    baf_plot, cn_plot, datatype_summary, maf, html_file, out_dir, sample_id
):

    cmd = ['run_report.sh', html_file, sample_id, mutations_per_cell, summary, 
        snvs_high_impact, snvs_all, trinuc, snv_adjacent_distance, snv_genome_count, 
        snv_cell_counts, snv_alt_counts, rearranegementtype_distribution_destruct_unfiltered, 
        chromosome_types_destruct_unfiltered,rearranegementtype_distribution_destruct_filtered, 
        chromosome_types_destruct_filtered,rearranegementtype_distribution_lumpy_unfiltered, 
        chromosome_types_lumpy_unfiltered,
        baf_plot, cn_plot, datatype_summary, maf]
    pypeliner.commandline.execute(*cmd)

    
def create_mutation_report(pseudobulk_group, merged_maf, high_impact_maf, high_impact_snvs, report_html):

    cmd = [ "run_mutationreport.sh", report_html,
        pseudobulk_group,
        high_impact_snvs,
        merged_maf,
        high_impact_maf]
    pypeliner.commandline.execute(*cmd)
