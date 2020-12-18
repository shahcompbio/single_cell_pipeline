import pypeliner
import pandas as pd
from single_cell.utils.csvutils import concatenate_csv, write_dataframe_to_csv_and_yaml
import matplotlib.pyplot as plt
import pandas as pd
import pypeliner
from single_cell.utils import helpers
from classifycopynumber import parsers, transformations
import os



def merge_segmental_cn(segmental_cn, concats):
    files = [pd.read_csv(f, sep="\t") for f in list(segmental_cn.values())]
    segmental_cn_combined = pd.concat(files)
    segmental_cn_combined.to_csv(concats, sep="\t", index=False)


def generate_segmental_copynumber(hmmcopy_files, segmental_cn, sample):


    cn, ploidy = parsers.read_hmmcopy_files(list(hmmcopy_files.values()), filter_normal=False, group_label_col='cell_id')

    cn["sample"] = sample
    transformations.generate_segmental_cn(segmental_cn, cn, ploidy, cn_col="copy", length_col="width")


def merge_cna_tables(tables, output):
    concatenate_csv(list(tables.values()), output)


def make_cna_table(amps, dels, cna_table):
    cna = {d1: [a,d] for (d1, a), (d2, d) in zip(amps.items(), dels.items()) }
    DATA=[]
    for label, files in cna.items():
        for f in files:
        
            data = pd.read_csv(f, usecols=["cn_type", "gene_name", "pass_filter"])
            data=data[data.pass_filter == True]
            data=data.rename(columns={"cn_type": "CN", "gene_name":"Gene"})
            data["Sample_Name"] = label
            # data["library"] = label
            # data["sample"] = sample
            data = data[["Gene", "Sample_Name", "CN"]]
            DATA.append(data)

    output = pd.concat(DATA)
    write_dataframe_to_csv_and_yaml(output, cna_table, output.dtypes, write_header=True)


def _write_maf(m, merged_maf, write_header):
    maf = pd.read_csv(m, sep="\t", chunksize=10e6)
    for chunk in maf:
        chunk.to_csv(merged_maf, sep="\t", index=False, header=write_header, mode='a')
        write_header=False   


def merge_mafs(germline_mafs, somatic_mafs, merged_maf):

    write_header=True

    for label, m in germline_mafs.items():
        _write_maf(m, merged_maf, write_header)
        write_header=False

    for label, m in somatic_mafs.items():
        _write_maf(m, merged_maf, write_header)
        write_header=False


def classify_hmmcopy(sample_label, hmmcopy_files, gtf, output_dir, amps, dels, docker_image=None):
    # gtf="/juno/work/shah/users/grewald/TWINS_NEW_DATA/WGS_REFERENCE/databases/Homo_sapiens.GRCh37.73.gtf"

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
    helpers.makedirs(tmpspace)

    cmd = [
        "/juno/work/shah/abramsd/oncokb-annotator/MafAnnotator.py", "-i", maf, "-o", annotated_maf, "-b", api_key
    ]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

def filter_maf(annotated_maf, filtered_maf, write_header=True):
    oncogenic_annotations = ["Oncogenic", "Likely Oncogenic", "Predicted Oncogenic"]
    maf = pd.read_csv(annotated_maf, sep="\t", chunksize=10e6)
    for chunk in maf:
        chunk = chunk[chunk.oncogenic.isin(oncogenic_annotations)]
        chunk.to_csv(filtered_maf, sep="\t", index=False, header=write_header, mode='a')

        write_header=False

def annotate_germline_somatic(filtered_maf, annotated_maf, is_germline):
    maf = pd.read_csv(filtered_maf, sep="\t")
    maf["is_germline"] = [is_germline] * len(maf)
    maf.to_csv(annotated_maf, sep="\t", index=False)


def label_germline_somatic(row):
    if row.is_germline == True:
        return row.Variant_Classification + "_" + "germline"
    return row.Variant_Classification + "_" + "somatic"
  

def prepare_maf_for_maftools(cohort_label, filtered_maf, prepared_maf, non_synonymous_labels, vcNames):
    '''
    filter on non synonymous labels
    add germline/somatic annotate to variant classification
    add group/patient labels
    write out vcNames
    '''
    maf = pd.read_csv(filtered_maf, sep="\t")
    maf = maf[maf.Variant_Classification.isin(non_synonymous_labels)]
    maf["Variant_Classification"] = maf.apply(lambda row: label_germline_somatic(row), axis=1 )
    nonsynclasses = pd.DataFrame({"Variant_Classification":maf.Variant_Classification.unique().tolist()})
    nonsynclasses.to_csv(vcNames, index=False)
    maf.to_csv(prepared_maf, sep="\t", index=False)


def make_oncoplot(prepped_maf, cna_table, oncoplot, vcNames, docker_image=None):

    plots_cmd = [
        "oncoplot.R", prepped_maf, vcNames, cna_table, oncoplot
    ]

    pypeliner.commandline.execute(*plots_cmd, docker_image=docker_image)

def format_cna_table(cna, out):
    data = pd.read_csv(cna)
    data.to_csv(out, sep="\t", index=False)
