import os
import pandas as pd
import pypeliner
from classifycopynumber import parsers, transformations
from single_cell.utils import helpers
import mafannotator.MafAnnotator as ma
import numpy as np
import shutil
import yaml
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
rmarkdown = importr("rmarkdown")


def merge_segmental_cn(segmental_cn, combined):
    '''
    merge copynumber in segmental format
    Parameters
    ----------
    segmental_cn : dictionary of sample: segmental cn
    combined : merged segmental cn data
    Returns
    -------
    '''
    files = [pd.read_csv(f, sep="\t") for f in list(segmental_cn.values())]
    segmental_cn_combined = pd.concat(files)
    segmental_cn_combined.to_csv(combined, sep="\t", index=False)

def generate_segmental_cn(filename, aggregated_cn_data, ploidy,  cn_col="copy", length_col="length"):

    aggregated_cn_data['ploidy'] = ploidy 
    aggregated_cn_data['seg.mean'] = np.log2(aggregated_cn_data[cn_col] / aggregated_cn_data['ploidy'])
    aggregated_cn_data['num.mark'] = (aggregated_cn_data[length_col] / 500000).astype(int)
    aggregated_cn_data = aggregated_cn_data.rename(columns={'sample': 'ID', 'chromosome': 'chrom', 'start': 'loc.start', 'end': 'loc.end'})
    aggregated_cn_data = aggregated_cn_data[['ID', 'chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean']]
    aggregated_cn_data['seg.mean'] = aggregated_cn_data['seg.mean'].fillna(np.exp(-8))
    aggregated_cn_data.loc[aggregated_cn_data['seg.mean'] == np.NINF, 'seg.mean'] = np.exp(-8)
    aggregated_cn_data = transformations._correct_seg_bin_ends(aggregated_cn_data)
    aggregated_cn_data.to_csv(filename, index=None, sep='\t')

def generate_segmental_copynumber(hmmcopy_files, segmental_cn, sample):
    '''
    transform copynumber  data from classify_copynumber
    package to segmental format
    Parameters
    ----------
    hmmcopy_files : dictionary of library: copynumber data
    segmental_cn : output segmental cn file, 1 for sample
    sample : sample label
    Returns
    -------
    '''
    cn, ploidy = parsers.read_hmmcopy_files(
        list(hmmcopy_files.values()),
        filter_normal=False, group_label_col='cell_id'
    )

    cn["sample"] = sample
    generate_segmental_cn(
        segmental_cn, cn, ploidy, cn_col="copy",
        length_col="width"
    )

def merge_vcfs_mafs(museq, strelka_snv, strelka_indel, merged, label):
    mafs = {"museq": museq, "strelka_snv": strelka_snv, "strelka_indel": strelka_indel}
    merge_mafs(mafs, merged, forced_label=label)


def merge_mafs(mafs, merged_maf, forced_label=None):
    """Write maf m to path merged_maf with label label (append).
    Args:
        mafs ([dict]): [sample: maf path]
        merged_maf ([str]): [output path]
    Yields:
        [type]: [description]
    """
    def _read_maf_chunk(maf_file, label):
        skiprows=0
        if open(maf_file).readline().startswith("#version"):
            skiprows=1
        maf = pd.read_csv(maf_file, sep="\t", dtype='str', chunksize=1e7, skiprows=skiprows)
        for maf_chunk in maf:
            if forced_label != None:
                maf_chunk["Tumor_Sample_Barcode"] = forced_label
            else:
                maf_chunk["Tumor_Sample_Barcode"] = label
            yield maf_chunk

    if os.path.exists(merged_maf):
        os.remove(merged_maf)

    for label, maf_file in mafs.items():
        for i, chunk in enumerate(_read_maf_chunk(maf_file, label)):
            header = True if i == 0 else False
            chunk.to_csv(
                merged_maf, sep="\t", index=False, header=header,
                mode='a', na_rep=""
            )


def generate_gistic_outputs(gistic_data, hdel_data, cbio_table):
    """Transform copynumber data to gistic format.
    Args:
        gistic_data ([str]): [input amp dataa]
        hdel_data ([str]): [input del ddata]
        cbio_table ([str]): [path to output table]
    """
    gistic_data['gistic_value'] = 2
    gistic_data.loc[gistic_data['log_change'] < 1, 'gistic_value'] = 1
    gistic_data.loc[gistic_data['log_change'] < 0.5, 'gistic_value'] = 0
    gistic_data.loc[gistic_data['log_change'] < -0.5, 'gistic_value'] = -1

    # Merge hdels
    hdel_data['is_hdel'] = 1
    gistic_data = gistic_data.merge(
        hdel_data[['Hugo_Symbol', 'sample', 'is_hdel']], how='left'
    )
    gistic_data['is_hdel'] = gistic_data['is_hdel'].fillna(0).astype(int)
    gistic_data.loc[gistic_data['is_hdel'] == 1, 'gistic_value'] = -2

    # Gistic_data generation
    gistic_data = gistic_data[['Hugo_Symbol', 'sample', 'gistic_value']]
    gistic_data = gistic_data.drop_duplicates()
    gistic_matrix = gistic_data.set_index(
        ['Hugo_Symbol', 'sample'])['gistic_value'].unstack()
    gistic_matrix.reset_index(inplace=True)

    gistic_matrix.to_csv(cbio_table, sep="\t", index=False)


def make_cbio_cna_table(cn_change_filename, cbio_table):
    gistic_data = pd.read_csv(cn_change_filename, sep="\t", usecols=["gene_name", "sample", "gistic_value"])
    gistic_data = gistic_data.rename(columns={"gene_name": "Hugo_Symbol"})
    gistic_data = gistic_data.drop_duplicates()
    gistic_data = gistic_data.astype({"gistic_value": "Int64"})

    gistic_matrix = gistic_data.set_index(['Hugo_Symbol', 'sample'])['gistic_value'].unstack()
    gistic_matrix.reset_index(inplace=True)
    gistic_matrix.to_csv(cbio_table, sep="\t", index=False, na_rep="NA")


def make_maftools_cna_table(cn_change_filename, maftools_table):
    cn_change = pd.read_csv(cn_change_filename, sep="\t",
                            usecols=["gene_name", "sample", "is_hdel", "is_loh", "is_hlamp"])

    cn_change = cn_change[cn_change['is_hdel'] | cn_change['is_loh'] | cn_change['is_hlamp']]

    cn_change["is_loh"] = cn_change.is_loh.replace({False: "", True: "loh"})
    cn_change["is_hdel"] = cn_change.is_hdel.replace({False: "", True: "hdel"})
    cn_change["is_hlamp"] = cn_change.is_hlamp.replace({False: "", True: "hlamp"})
    cn_change["CN"] = cn_change[["is_hdel", "is_loh", "is_hlamp"]].apply(
        lambda row: "".join(map(str, row)) ,axis=1
    )

    cn_change = cn_change.rename(columns={"gene_name": "Gene", "cn_type": "CN", "sample": "Sample_Name"})
    cn_change = cn_change[["Gene", "Sample_Name", "CN"]]
    cn_change.to_csv(maftools_table, index=False, sep="\t")


def merge_cna_tables(tables, output):
    '''
    merge copynumber tables from classify_copynumber
    Parameters
    ----------
    tables : sample: path of cn data from classify copynumber
    output: output
    Returns
    -------
    '''
    number = 0
    for label, cna in tables.items():

        data = pd.read_csv(cna)
        data["sample"] = label
        if number == 0:
            header = True
        else:
            header = False
        number += 1
        data.to_csv(output, index=False, mode='a', header=header, sep="\t")


def _get_hmmcopy_sample_ids(files):
    samples = []
    for library, f in files.items():
        mf = yaml.load(open(os.path.join(os.path.dirname(f), "metadata.yaml")))
        ids = mf["meta"]["cell_ids"]
        samples += [id.split("-")[0] for id in ids]
    return pd.Series(samples).unique().tolist()


def classify_hmmcopy(
    hmmcopy_files, gtf, output_dir,
    amps, dels, docker_image=None
):
    '''
    run classify_copynumber on hmmcopy data
    Parameters
    ----------
    sample_label: sample label for all libraries
    hmmcopy_files: library: path of hmmcopy data 
    gtf: gtf file
    output_dir: output directory
    amps: output path for amps
    dels: output path for dels
    Returns
    -------
    '''
    sample_ids = _get_hmmcopy_sample_ids(hmmcopy_files)
    files = list(hmmcopy_files.values())
    cmd = [
        "classifycopynumber", gtf, output_dir,
         amps, dels
    ]
    for f in files:
        cmd.extend(["--hmmcopy_csv_filenames", f])
    for s_id in sample_ids:
        cmd.extend(["--sample_ids", s_id])

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def annotate_maf_with_oncokb(
        maf, api_key, tmpspace, annotated_maf
):
    '''
    annotate maf with onco kb to get oncogenicity of variants
    Parameters
    ----------
    maf: input maf
    api_key: api key for oncokb api 
    tmpspace: tmp dir to run in 
    annotated_maf: outputed annotated dmaf
    Returns
    -------
    '''
    helpers.makedirs(tmpspace)
    ma.annotate(maf, annotated_maf, api_key)


def filter_maf(annotated_maf, filtered_maf, write_header=True):
    '''
    filter maf on annotated oncogenicity
    Parameters
    ----------
    annotated_maf: annotated maf
    filtered_maf: output filtered_maf
    write_header: bool, write first  header
    Returns
    -------
    '''
    oncogenic_annotations = [
        "Oncogenic", "Likely Oncogenic", "Predicted Oncogenic"
    ]
    maf = pd.read_csv(annotated_maf, sep="\t", dtype='str', chunksize=10e6)
    for chunk in maf:
        chunk = chunk[chunk.ONCOGENIC.isin(oncogenic_annotations)]
        chunk.to_csv(
            filtered_maf, sep="\t", index=False, header=write_header, mode='a'
        )
        write_header = False


def annotate_maf_file(filtered_maf, annotated_maf, annotations):
    maf = pd.read_csv(filtered_maf, sep="\t", dtype='str')

    for k, v in annotations.items():
        assert isinstance(v, int) or isinstance(v, float) or isinstance(v, str)
        maf[k] = v

    maf.to_csv(annotated_maf, sep="\t", index=False)


def label_germline_somatic(row):
    '''
    add germline/somatic label to row  of maf
    Parameters
    ----------
    row: row of pandas df
    Returns
    -------
    '''
    if row.is_germline == True:
        return row.Variant_Classification + "_" + "germline"
    return row.Variant_Classification + "_" + "somatic"


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
        vep_ref['reference_dir'],
        # "--buffer_size",
        # "20"
    ]

    pypeliner.commandline.execute(*cmd)


def prepare_maf_for_maftools(
    germline, somatic, prepared_maf, non_synonymous_labels, vcNames
):
    '''
    format maf for intake in maftools
    Parameters
    ----------
    filtered_maf: label of filtered_maf
    prepared_maf: new formatted maf
    non_synonymous_labels: set of non_synonymous variant types to include
    vcNames: set of non_synonymous variant types existing in maf, passed to maftools
    Returns
    -------
    '''
    germline = pd.read_csv(germline, sep="\t")
    germline = germline[
        germline.Variant_Classification.isin(non_synonymous_labels)
    ]
    germline["Variant_Classification"] = germline.Variant_Classification.apply(
        lambda vc: vc + "_germline"
    )

    somatic = pd.read_csv(somatic, sep="\t")
    somatic = somatic[
        somatic.Variant_Classification.isin(non_synonymous_labels)
    ]
    somatic["Variant_Classification"] = somatic.Variant_Classification.apply(
        lambda vc: vc + "_somatic"
    )

    combined = pd.concat([germline, somatic])

    nonsynclasses = pd.DataFrame(
        {"Variant_Classification": combined.Variant_Classification.unique().tolist()}
    )
    nonsynclasses.to_csv(vcNames, index=False)

    combined.to_csv(prepared_maf, sep="\t", index=False)


def make_oncoplot(
    prepped_maf, cna_table, oncoplot, vcNames, docker_image=None
):
    '''
    run R script to make oncoplot with prepped, annotatedd maf
    Parameters
    ----------
    prepped_maf: formatted, annotated maf
    cna_table: cna data formatted for maftools
    oncoplot: path for oncplot
    vcNames: set of non_synonymous variant types existing in maf, passed to maftools
    Returns
    -------
    '''
    script_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'scripts','oncoplot.R')

    plots_cmd = [
        "Rscript", script_path, prepped_maf, vcNames, cna_table, oncoplot
    ]
    pypeliner.commandline.execute(*plots_cmd, docker_image=docker_image)


def create_report(
        cohort, oncoplot, report
):
    rmd_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        'scripts', 'report.Rmd'
    )

    parameters = robjects.r.list(cohort=cohort, oncoplot=oncoplot)

    rmarkdown.render(rmd_script, output_file=report,
        output_options=robjects.r.list(self_contained=True),
        params=parameters
    )
