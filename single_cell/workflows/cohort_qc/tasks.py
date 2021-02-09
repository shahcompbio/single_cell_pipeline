import pypeliner
import pandas as pd
from single_cell.utils.csvutils import concatenate_csv, write_dataframe_to_csv_and_yaml
import matplotlib.pyplot as plt
import pandas as pd
import pypeliner
from single_cell.utils import helpers
from classifycopynumber import parsers, transformations
import os


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
    cn, ploidy = parsers.read_hmmcopy_files(list(hmmcopy_files.values()), 
        filter_normal=False, group_label_col='cell_id')

    cn["sample"] = sample
    transformations.generate_segmental_cn(segmental_cn, cn, ploidy, cn_col="copy", 
        length_col="width")


def _write_maf(m, label, merged_maf, write_header):
    '''
    write maf m to path merged_maf with label label (append)
    Parameters
    ----------
    m : maf path
    labe: maf label
    merged_maf : write path
    write_header: bool, write header
    Returns
    -------
    '''
    maf = pd.read_csv(m, sep="\t", dtype='str', chunksize=10e6)
    for chunk in maf:
        chunk["Tumor_Sample_Barcode"] = label
        chunk= chunk.astype({"t_ref_count":"Int64", "t_alt_count":"Int64", 
            "n_ref_count":"Int64", "n_alt_count":"Int64"})

        chunk.to_csv(merged_maf, sep="\t", index=False, header=write_header, mode='a', na_rep="")
        write_header=False     


def merge_mafs(germline, somatic_mafs, merged_maf):
    '''
    write maf m to path merged_maf with label label (append)
    Parameters
    ----------
    m : maf path
    labe: maf label
    merged_maf : write path
    write_header: bool, write header
    Returns
    -------
    '''
    write_header=True

    for label, m in somatic_mafs.items():
        _write_maf(m, label, merged_maf, write_header)
        write_header=False
    for label, m in germline.items():
        _write_maf(m, label, merged_maf, write_header)
        write_header=False


def generate_gistic_outputs(gistic_data, hdel_data, cbio_table):
    '''
    transform copynumber data to gistic format
    Parameters
    ----------
    gistic_data : amps
    hdel_data: dels
    cbio_table: output
    Returns
    -------
    '''
    gistic_data['gistic_value'] = 2
    gistic_data.loc[gistic_data['log_change'] < 1, 'gistic_value'] = 1
    gistic_data.loc[gistic_data['log_change'] < 0.5, 'gistic_value'] = 0
    gistic_data.loc[gistic_data['log_change'] < -0.5, 'gistic_value'] = -1

    # Merge hdels
    hdel_data['is_hdel'] = 1
    gistic_data = gistic_data.merge(hdel_data[['Hugo_Symbol', 'sample', 'is_hdel']], how='left')
    gistic_data['is_hdel'] = gistic_data['is_hdel'].fillna(0).astype(int)
    gistic_data.loc[gistic_data['is_hdel'] == 1, 'gistic_value'] = -2

    # Gistic_data generation
    gistic_data = gistic_data[['Hugo_Symbol', 'sample', 'gistic_value']]
    gistic_data = gistic_data.drop_duplicates()
    gistic_matrix = gistic_data.set_index(['Hugo_Symbol', 'sample'])['gistic_value'].unstack()
    gistic_matrix.reset_index(inplace=True)

    gistic_matrix.to_csv(cbio_table, sep="\t", index=False)


def make_cbio_cna_table(amps, dels, cbio_table):
    '''
    transform amps, dels to cbio-readable format
    Parameters
    ----------
    amps : amps
    dels: dels
    cbio_table: output
    Returns
    -------
    '''
    amps = pd.read_csv(amps,  sep="\t", usecols=["gene_name", "log_change", "sample"])
    amps = amps.rename(columns={"gene_name":"Hugo_Symbol"})

    dels = pd.read_csv(dels,  sep="\t", usecols=["gene_name", "sample"])
    dels = dels.rename(columns={"gene_name":"Hugo_Symbol"})

    generate_gistic_outputs(amps, dels, cbio_table)


def make_maftools_cna_table(amps, dels, maftools_table):
    '''
    transform amps, dels to maftools-readable format
    Parameters
    ----------
    amps : amps
    dels: dels
    maftools_table: output
    Returns
    -------
    '''
    amps = pd.read_csv(amps, sep="\t", usecols=["gene_name", "sample", "cn_type", "pass_filter"])
    amps = amps.rename(columns={"gene_name":"Gene", "cn_type":"CN", "sample": "Sample_name"})
    amps=amps[amps.pass_filter == True]

    dels = pd.read_csv(dels,  sep="\t", usecols=["gene_name", "sample", "cn_type", "pass_filter"])
    dels = dels.rename(columns={"gene_name":"Gene", "cn_type":"CN", "sample": "Sample_name"})
    dels=dels[dels.pass_filter == True]

    out = pd.concat([amps, dels])
    out = out[["Gene","Sample_name","CN"]]
    out.to_csv(maftools_table, index=False, sep="\t")


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
    number=0
    for label, cna in tables.items():

        data = pd.read_csv(cna)
        data["sample"] = label
        if number==0:
            header=True
        else:
            header=False
        number+=1
        data.to_csv(output, index=False, mode='a', header=header, sep="\t")


def classify_hmmcopy(sample_label, hmmcopy_files, gtf, output_dir, amps, dels, docker_image=None):
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
    files = list(hmmcopy_files.values())

    cmd = [
        "classifycopynumber", gtf, output_dir, sample_label, amps, dels, "--plot", False
    ]
    for f in files:
        cmd.extend(["--hmmcopy_csv_filenames", f])

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def annotate_maf_with_oncokb(

        maf, api_key, tmpspace, annotated_maf, docker_image=None
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

    cmd = [
        "MafAnnotator.py", "-i", maf, "-o", annotated_maf, "-b", api_key
    ]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


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
    oncogenic_annotations = ["Oncogenic", "Likely Oncogenic", "Predicted Oncogenic"]
    maf = pd.read_csv(annotated_maf, sep="\t", dtype='str', chunksize=10e6)
    for chunk in maf:
        chunk = chunk[chunk.oncogenic.isin(oncogenic_annotations)]
        chunk.to_csv(filtered_maf, sep="\t", index=False, header=write_header, mode='a')

        write_header=False


def annotate_germline_somatic(filtered_maf, annotated_maf, is_germline):
    '''
    add germline/somatic annotated to filtered_maf
    Parameters
    ----------
    filtered_maf: filtered_maf maf
    annotated_maf: output annotated_maf
    is_germline: str, label to add to maf, should  be in {germline, somatic}
    Returns
    -------
    '''
    maf = pd.read_csv(filtered_maf, sep="\t", dtype='str')
    maf["is_germline"] = [is_germline] * len(maf)
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
  

def prepare_maf_for_maftools(cohort_label, filtered_maf, prepared_maf, non_synonymous_labels, vcNames):
    '''
    format maf for intake in maftools
    Parameters
    ----------
    cohort_label: label of cohort
    filtered_maf: label of filtered_maf
    prepared_maf: new formatted maf
    non_synonymous_labels: set of non_synonymous variant types to include
    vcNames: set of non_synonymous variant types existing in maf, passed to maftools

    Returns
    -------
    '''
    maf = pd.read_csv(filtered_maf, sep="\t", dtype='str')
    maf = maf[maf.Variant_Classification.isin(non_synonymous_labels)]
    maf["Variant_Classification"] = maf.apply(lambda row: label_germline_somatic(row), axis=1 )
    nonsynclasses = pd.DataFrame({"Variant_Classification":maf.Variant_Classification.unique().tolist()})
    nonsynclasses.to_csv(vcNames, index=False)
    maf.to_csv(prepared_maf, sep="\t", index=False)


def make_oncoplot(prepped_maf, cna_table, oncoplot, vcNames, docker_image=None):
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
    plots_cmd = [
        "oncoplot.R", prepped_maf, vcNames, cna_table, oncoplot
    ]
    pypeliner.commandline.execute(*plots_cmd, docker_image=docker_image)
