import os
import pandas as pd
import pypeliner
from classifycopynumber import parsers, transformations
from single_cell.utils import helpers


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


def generate_segmental_copynumber(hmmcopy_files, sample_ids, segmental_cn, sample_label):
    '''
    transform copynumber  data from classify_copynumber
    package to segmental format
    Parameters
    ----------
    hmmcopy_files : dictionary of library: copynumber data
    sample_ids : sample_ids for filtering hmmcopy data
    segmental_cn : output segmental cn file, 1 for sample
    sample_label : sample label for merged data
    Returns
    -------
    '''
    sample_ids = list(sample_ids.values())
    cn, ploidy = parsers.read_hmmcopy_files(
        list(hmmcopy_files.values()),
        filter_normal=False, group_label_col='cell_id',
        sample_ids=sample_ids,
    )

    cn["sample"] = sample_label
    transformations.generate_segmental_cn(
        segmental_cn, cn, ploidy, cn_col="copy",
        length_col="width"
    )


def merge_mafs(mafs, merged_maf):
    """Write maf m to path merged_maf with label label (append).

    Args:
        mafs ([dict]): [sample: maf path]
        merged_maf ([str]): [output path]

    Yields:
        [type]: [description]
    """
    def _read_maf_chunk(maf_file, label):
        maf = pd.read_csv(maf_file, sep="\t", dtype='str', chunksize=1e7)
        for maf_chunk in maf:
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
    amps = pd.read_csv(
        amps, sep="\t", usecols=["gene_name", "log_change", "sample"]
    )
    amps = amps.rename(columns={"gene_name": "Hugo_Symbol"})

    dels = pd.read_csv(dels, sep="\t", usecols=["gene_name", "sample"])
    dels = dels.rename(columns={"gene_name": "Hugo_Symbol"})

    generate_gistic_outputs(amps, dels, cbio_table)


def make_maftools_cna_table(amps, dels, maftools_table):
    """Transform amps, dels to maftools-readable format/

    Args:
        amps ([dict]): [amps]
        dels ([dict]): [dels]
        maftools_table ([str]): [output path]
    """
    amps = pd.read_csv(
        amps, sep="\t", usecols=["gene_name", "sample",
             "cn_type", "pass_filter"]
    )
    amps = amps.rename(
        columns={"gene_name": "Gene", "cn_type": "CN",
            "sample": "Sample_name"}
        )
    amps = amps[amps.pass_filter == True]

    dels = pd.read_csv(
        dels,  sep="\t", usecols=["gene_name", "sample",
            "cn_type", "pass_filter"]
    )
    dels = dels.rename(
        columns={"gene_name": "Gene", "cn_type": "CN", "sample": "Sample_name"}
    )
    dels = dels[dels.pass_filter == True]

    out = pd.concat([amps, dels])
    out = out[["Gene", "Sample_name", "CN"]]
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


def classify_hmmcopy(
    hmmcopy_files, sample_ids, gtf, output_dir, 
    amps, dels
):
    '''
    run classify_copynumber on hmmcopy data
    Parameters
    ----------
    hmmcopy_files: library: path of hmmcopy data 
    sample_ids: sample ids to filter library level data
    gtf: gtf file
    output_dir: output directory
    amps: output path for amps
    dels: output path for dels
    Returns
    -------
    '''
    files = list(hmmcopy_files.values())
    sample_ids = list(sample_ids.values())
    cmd = [
        "classifycopynumber", gtf, output_dir,
        amps, dels, "--plot", False
    ]
    for f in files:
        cmd.extend(["--hmmcopy_csv_filenames", f])
    for s in sample_ids:
        cmd.extend(["--sample_ids", s])

    pypeliner.commandline.execute(*cmd)


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

    cmd = [
         "MafAnnotator.py", "-i", maf, "-o", annotated_maf, "-b", api_key
    ]

    pypeliner.commandline.execute(*cmd)


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
    prepped_maf, cna_table, oncoplot, vcNames
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
    plots_cmd = [
        "oncoplot.R", prepped_maf, vcNames, cna_table, oncoplot
    ]
    pypeliner.commandline.execute(*plots_cmd)
