import os

import numpy as np
import pandas as pd
import pypeliner
from single_cell.utils import csvutils
from single_cell.utils import helpers


def parse_vcf(vcf, csv, return_pandas=False):
    data = []
    cols = []
    with open(vcf) as fp:
        for cnt, line in enumerate(fp):
            if line[0] != "#":
                data.append(line.strip().split("\t"))
            if line[0] == "#" and line[1] != "#":
                cols = line.strip().split("\t")
                if cols[0] == "#CHROM": cols[0] = "CHROM"
    vcf = pd.DataFrame(data, columns=cols)
    if return_pandas:
        return vcf
    if not return_pandas:
        csv = pd.read_csv(csv, delimiter="\t")
        csv.to_csv(csv, index=False, encoding="utf-8-sig")


class VarcallLoader:
    def __init__(self, csv, caller):

        self.caller = caller

        assert caller in ['lumpy', 'destruct']
        delim = ',' if caller == 'lumpy' else '\t'

        self.svtyper_set = [
            "CHROM", "POS", "CHROM2", "POS2"
        ]

        self.lumpy_set = [
            "chrom1", "start1", "chrom2",
            "start2"
        ]

        self.destruct_set = [
            "chromosome_1", "position_1",
            "chromosome_2",
            "position_2"
        ]

        self.lumpy_to_svtype = {
            "INTERCHROM": "BND",
            "DUPLICATION": "DUP",
            "DELETION": "DEL",
            "INSERTION": "INS",
            "INVERSION": "INV"
        }

        self.destruct_to_svtype = {
            "translocation": "BND",
            "duplication": "DUP",
            "deletion": "DEL",
            "insertion": "INS",
            "inversion": "INV"
        }

        self.data = pd.read_csv(csv, delim)

        self.translator = dict(zip(self.svtyper_set, self.destruct_set))

        if caller == "lumpy":
            self.translator = dict(zip(self.svtyper_set, self.lumpy_set))

    def __getitem__(self, index):
        if index == "STRAND" and self.caller == "lumpy":
            return self.data["strands"]

        if index == "STRAND" and self.caller == "destruct":
            return self.data["strand_1"].combine(self.data["strand_2"],
                                                 lambda s1, s2: str(str(s1) + str(s2)))

        if index == "TYPE" and self.caller == "lumpy":
            return self.data["type"].apply(lambda type: self.lumpy_to_svtype[type])

        if index == "TYPE" and self.caller == "destruct":
            return self.data["type"].apply(lambda type: self.destruct_to_svtype[type])

        else:
            return self.data[self.translator[index]]


def make_alt(strand, pos, chrom, ref):
    """
    create an alt section for an output vcf
    :param strand: string representing the strand senses
    :type strand: str
    :param pos: pos of alt
    :type pos: int
    :param chrom: chrom of alt
    :type chrom: str
    :param ref: reference seq
    :type ref: str
    :return: vcf format alt
    :rtype: str
    """
    insert = str(chrom) + ":" + str(pos)
    ref = str(ref)
    if strand == "--":
        return "[" + insert + "[" + ref
    elif strand == "++":
        return ref + "]" + insert + "]"
    elif strand == "-+":
        return "]" + insert + "]" + ref
    elif strand == "+-":
        return ref + "[" + insert + "["
    else:
        return 0


def expand_info_section(svtype, loaderobj, index):
    """
    creates the info section for
    the outputvcf from input destruct
    csv
    :param svtype: type of structural variation; see dict "to_svtype" in csv2vcf() for options
    :type svtype: str
    :param loaderobj:
    :type loaderobj:
    :param index: current_index
    :type index: str
    :param variant_caller: caller (destruct/lumpy)
    :type variant_caller: str
    :return: info sectiopn
    :rtype: dict
    """
    end = loaderobj["POS2"][index]

    start_ci = "0,200"
    end_ci = "0,200"

    info = {'SVTYPE': svtype,
            'CIPOS': start_ci,
            'CIPOS95': start_ci,
            'CIEND': end_ci,
            'CIEND95': end_ci}

    if svtype is "BND":
        info["MATEID"] = str(index) + "_2"

    else:
        info["END"] = end
    return info


def info_tostr(info):
    if type(info) == dict:
        return ";".join(("{}={}".format(*i) for i in info.items()))
    else:
        return info


def add_row(data, new_row, index):
    before = data.iloc[:index]
    after = data.iloc[index:]

    return pd.concat([before, new_row, after]).reset_index(drop=True)


def add_bnd_mates(new_data, loaderobj, num_records, variant_caller):
    # loop through everything
    row_num = 0
    old_data_index = 0

    while row_num < num_records:
        if new_data["INFO"][row_num]["SVTYPE"] == "BND":
            chrom = loaderobj["CHROM2"][old_data_index]
            pos = loaderobj["POS2"][old_data_index]
            id = str(old_data_index) + "_2"
            ref = new_data["REF"][row_num]
            strand = loaderobj["STRAND"][old_data_index]
            strand = strand[1] + strand[0]
            alt = make_alt(strand, pos, chrom, ref)
            qual = "."
            vcf_filter = new_data["FILTER"][row_num]
            info = expand_info_section("BND", loaderobj, old_data_index)
            info["MATEID"] = str(old_data_index) + "_1"

            bnd_mate = pd.DataFrame({'#CHROM': chrom, 'POS': pos, 'ID': id,
                                     'REF': ref, 'ALT': alt, 'QUAL': qual,
                                     'FILTER': vcf_filter,
                                     'INFO': info_tostr(info)},
                                    index=[row_num],
                                    columns=['#CHROM', 'POS', 'ID',
                                             'REF', 'ALT', 'QUAL',
                                             'FILTER',
                                             'INFO'])

            new_data = add_row(new_data, bnd_mate, row_num + 1)

            num_records += 1
            row_num += 1
        row_num += 1
        old_data_index += 1

    return new_data


def varcalls_to_svtyper_input(input, vcf, tempdir, caller):
    """
    convert input variant call file csv
    [either lumpy or destruct]
    to vcf that svtyper can take

    :param input: input variant call CSV
    :type input: string
    :param vcf: output vcf
    :type vcf: string
    :param tempdir: temp dir
    :type tempdir: string
    :param caller: variant caller, either lumpy or destruct
    :type caller: string
    """

    helpers.makedirs(tempdir)

    csv = VarcallLoader(input, caller)
    num_records = len(csv.data.index)

    chrom = csv["CHROM"]

    pos = csv["POS"]

    ref = ["N"] * num_records

    qual = ["."] * num_records

    vcf_filter = ["Pass"] * num_records

    sv_type = csv["TYPE"]

    vcf_id = [None] * num_records
    alt = [None] * num_records
    info = [None] * num_records

    for i in range(num_records):
        if sv_type[i] == "BND":
            vcf_id[i] = str(i) + "_1"
            alt[i] = make_alt(
                csv["STRAND"][i],
                pos[i],
                chrom[i],
                ref[i]
            )
        else:
            vcf_id[i] = str(i)
            alt[i] = "<" + str(sv_type[i]) + ">"

        info[i] = expand_info_section(
            sv_type[i], csv, i
        )

        new_data = pd.DataFrame(
            {'#CHROM': chrom, 'POS': pos, "ID": vcf_id,
             "REF": ref, 'ALT': alt, 'QUAL': qual,
             'FILTER': vcf_filter, 'INFO': info
             },
            columns=[
                '#CHROM', 'POS', 'ID', 'REF',
                'ALT', 'QUAL', 'FILTER', 'INFO'
            ]
        )

    new_data = add_bnd_mates(new_data, csv, num_records, caller)

    new_data["INFO"] = [info_tostr(info) for info in new_data["INFO"]]

    new_data.to_csv(vcf, sep="\t", index=False)


def extract_svtyper_info(df):
    num_records = len(df.index)
    data = []
    headers = None
    for i in range(num_records):
        headers = df.iloc[i][-2].split(":")
        data.append(df.iloc[i][-1].split(":"))
    out = pd.DataFrame(data, columns=headers)
    return out


def genotype(
        input_bam, reference, input_vcf,
        output_vcf, output_csv, tempdir,
        cell_id, docker_image=None
):
    """
    calls svtyper-sso on input
    bam and vcf to perform genotyping.
    :param input_bam:
    :type input_bam:
    :param reference:
    :type reference:
    :param input_vcf:
    :type input_vcf:
    :param output_vcf:
    :type output_vcf:
    :param output_csv:
    :type output_csv:
    :param tempdir:
    :type tempdir:
    :param docker_image:
    :type docker_image:
    :return:
    :rtype:
    """
    helpers.makedirs(tempdir)

    cmd = ['svtyper-sso',
           '--input_vcf', input_vcf,
           '--bam', input_bam,
           '--ref_fasta', reference,
           '-o', output_vcf]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    base_data = parse_vcf(output_vcf, None, return_pandas=True)

    svtype_annotations = extract_svtyper_info(base_data)

    base_data = base_data.iloc[:, :-2]  # assumes svtyper info in last 2 cols

    output = pd.concat([base_data, svtype_annotations], axis=1)

    output['cell_id'] = cell_id

    csvutils.write_dataframe_to_csv_and_yaml(output, output_csv, {}, write_header=True)


def write_svtyper_annotation(annotation, inputcsv, outfile):
    """
    takes an inputcsv with an svtyper
    annotation and writes it out to outfile
    in matrix form.
    :param annotation: svtyper annotation to write
    :type annotation: str
    :param inputcsv: inputcsv containing svtyper annotation
    :type inputcsv: pd.DataFrame
    :param outfile: filename for write out file
    :type outfile: str
    """
    inputcsv = inputcsv[['CHROM', 'POS', annotation, 'cell_id']]

    inputcsv['coord'] = list(zip(inputcsv['CHROM'], inputcsv['POS']))

    inputcsv = inputcsv.pivot(index='cell_id', columns='coord', values=annotation)

    inputcsv = inputcsv.T

    inputcsv[['CHROM', 'POS']] = pd.DataFrame(inputcsv.index.tolist(), index=inputcsv.index)

    inputcsv.to_csv(outfile, index=False)


def write_svtyper_annotations(csv, output_paths, tempdir):
    """
    writes the annotations contained in the below
    annotations list to files, each to their own

    :param csv: csv file containg annotations as features
    :type csv:
    :param output_paths: output directories for annotation files
    :type output_paths:
    :param tempdir:
    :type tempdir:
    :return:
    :rtype:
    """
    helpers.makedirs(tempdir)

    annotations = ["AO", "AP", "AS",
                   "ASC", "DP", "GQ",
                   "QA", "QR", "RO",
                   "RP", "RS", "SQ",
                   "GL", "AB"]

    csv = pd.read_csv(csv, delimiter=",")

    for annotation in annotations:
        temp_output_path = os.path.join(tempdir, '{}.csv.gz'.format(annotation))

        write_svtyper_annotation(
            annotation,
            csv,
            temp_output_path
        )

        csvutils.rewrite_csv_file(temp_output_path, output_paths[annotation])


def merge_csvs(input_csvs, merged_csv):
    """
    merges input csv files
    into one csv
    """
    csvutils.concatenate_csv(input_csvs, merged_csv, write_header=True)
