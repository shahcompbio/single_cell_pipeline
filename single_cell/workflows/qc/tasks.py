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
    print(snv_files)
    
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
    cmd = ["Rscript", "/juno/work/shah/abramsd/CODE/single_cell_pipeline/single_cell/workflows/qc/scripts/mergesnvs.R", snv, filtsnv]
    pypeliner.commandline.execute(*cmd, docker_image="rocker/tidyverse")

def filter_maf_for_high_impact(maf, filtmaf):
    cmd = ["Rscript", "/juno/work/shah/abramsd/CODE/single_cell_pipeline/single_cell/workflows/qc/scripts/mergemafs.R", maf, filtmaf]
    pypeliner.commandline.execute(*cmd, docker_image="rocker/tidyverse")

def _get_indelvcfs(dir):
    indel_vcfs = []
    indel_vcfs += sorted(glob.glob(dir + "/results/variants/" + "*_strelka_indel.vcf.gz" ))
    indel_vcfs += sorted(glob.glob(dir + "/results/" + "*_strelka_indel.vcf.gz" ))
    indel_vcfs += sorted(glob.glob(dir + "/results/variant_calling/*/" + "*indel.vcf.gz" ))
    return indel_vcfs

#dummie outputs for testing
def dummie(files):
    for file in files:
        print(file)
        cc
        df = pd.DataFrame([1], [1])
        df.to_csv(file)
    return

def get_snv_all_csvs(dir):
    output = {}
    dir_len = len(dir.split("/"))
    filenames =  sorted(glob.glob(dir + "/*/*/*snv*all*csv", recursive=True))
    print(filenames)
    print(glob.glob(dir + "/*/*/*", recursive=True))
    print(dir)
    for filename in filenames:

        name_len = len(filename.split("/"))
        samplelabel = filename.split("/")[dir_len : name_len - 1]
        output[tuple(samplelabel)] = filename
        print(output)
    return output

def gunzip_file(infile, outfile):
    with gzip.open(infile, 'rb') as f_in:
        with open(outfile, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

# def merge_mafs(mafs, mergedmaf):

#     cmd = " awk '(NR == 2) || (FNR > 2)' {} > {}".format(*mafs, mergedmaf)
#     pypeliner.commandline.execute(cmd, docker_image="shub://rdmorin/cancer_docker_singularity:vcf2maf")

def vcf2maf(vcf_file, output_maf, tempdir, reference, vepdata):

    if vcf_file.endswith('.gz'):
        helpers.makedirs(tempdir)
        vcf_unzipped = os.path.join(tempdir, 'unzipped_vcf.vcf')
        gunzip_file(vcf_file, vcf_unzipped)
    else:
        vcf_unzipped = vcf_file

    cmd = ["vcf2maf.pl", "--input-vcf", vcf_unzipped, "--output-maf", output_maf,
           "--ref-fasta", reference, "--filter-vcf", "0", "--vep-path", "/home/abramsd/miniconda3/envs/r-environment/bin", 
           "--vep-data", vepdata]
    print(cmd)

    pypeliner.commandline.execute(*cmd, docker_image="shub://rdmorin/cancer_docker_singularity:vcf2maf")


def sample_level_report(mutations_per_cell, summary, 
                      snvs_high_impact, snvs_all, trinuc, snv_adjacent_distance, snv_genome_count, 
                      snv_cell_counts, snv_alt_counts, rearranegementtype_distribution, chromosome_types,
                      BAFplot, CNplot, datatype_summary, maf, html_file, out_dir, sample_id):


	rcode = "rmarkdown::render('/juno/work/shah/abramsd/CODE/single_cell_pipeline/single_cell/workflows/qc/scripts/report2.Rmd',output_file = '{}',  params=list( sample_id='{}', mutations_per_cell_png='{}', summary_csv='{}', snvs_high_impact_csv='{}', snvs_all_csv='{}', trinuc_csv='{}', snv_adjacent_distance_png='{}', snv_genome_count_png='{}', snv_cell_counts_png='{}', snv_alt_counts_png='{}', rearranegementtype_distribution_png='{}', chromosome_types_png='{}', BAFplot_png='{}', cn_plot_png='{}', datatype_summary_csv='{}', maf='{}'))".format(
					html_file, sample_id, mutations_per_cell, summary, 
					snvs_high_impact, snvs_all, trinuc, snv_adjacent_distance, snv_genome_count, 
					snv_cell_counts, snv_alt_counts, rearranegementtype_distribution, chromosome_types,
					BAFplot, CNplot, datatype_summary, maf)
	cmd = ["R", "-e",  rcode]
	pypeliner.commandline.execute(*cmd, docker_image="rocker/tidyverse")

    
def create_mutation_report(pseudobulk_group, merged_maf, high_impact_maf, high_impact_snvs, report_html):

    rcode = "rmarkdown::render('/juno/work/shah/abramsd/CODE/single_cell_pipeline/single_cell/workflows/qc/scripts/mutationreport.Rmd',output_file = '{}',  params=list(pseudobulk_group='{}', merged_filt_snvs='{}', merged_maf='{}',high_impact_maf='{}'))".format(
        report_html,
        pseudobulk_group,
        high_impact_snvs,
        merged_maf,
        high_impact_maf
    )
    cmd = ["R", "-e",  rcode]
    pypeliner.commandline.execute(*cmd, docker_image="rocker/tidyverse")
