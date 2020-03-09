import sys
from single_cell.utils import inpututils
import os

def load(yaml):
    """
    load a yaml
    :param yaml: yaml file
    :return: dict
    """
    return inpututils.load_yaml(yaml)

def could_be_file(f):
    """
    check if str f could be a file
    :param f: str
    :return: T/F
    """
    if not isinstance(f, str):
        return False
    return os.path.dirname(f) != "" \
           and os.path.basename(f).split(".")[-1] != "" \
           and os.path.basename(f) != ""


def flatten(d, sep="_"):
    """
    flatten dict d
    :param d: dict
    :param sep: sep for compound keys
    :return: ordereddict
    """
    import collections

    obj = collections.OrderedDict()

    def recurse(t, parent_key=""):

        if isinstance(t, list):
            for i in range(len(t)):
                recurse(t[i], parent_key + sep + str(i) if parent_key else str(i))
        elif isinstance(t, dict):
            for k, v in t.items():
                recurse(v, parent_key + sep + k if parent_key else k)
        else:
            obj[parent_key] = t

    recurse(d)

    return obj

def make_parsed_format_consistent(yaml):
    """
    make parse format consistent
    :param yaml: yaml
    :return:
    """
    reformmated_yaml = yaml

    if isinstance(yaml, tuple):
        reformmated_yaml = yaml[0]

    if isinstance(yaml, str):
        reformmated_yaml = {"file": yaml}

    return reformmated_yaml

def check_local_input_files(input_files):
    """
    check that all files exist
    :param input_files:
    :return:
    """
    for file in input_files:
        assert os.path.exists(file)
        assert os.access(file, os.R_OK) #readable

def check_local_output_file(output_files):
    """
    check output files can be written to
    :param output_files:
    :return:
    """
    for file in output_files:
        assert os.access(file, os.W_OK) #writeable

def validate_alignment_types(alignment_yaml, alignment_type_example):
    """
    make sure alignment input yaml has correct typing
    :param alignment_yaml:
    :param alignment_type_example:
    :return:
    """
    if alignment_yaml.keys() != alignment_type_example.keys():
        return False

    for key in alignment_type_example.keys():
        if type(alignment_type_example[key]) != type(alignment_yaml[key]):
            return False

    return True

def validate_yaml(yaml, type, local=True):
    """
    validate yaml
    :param yaml: yaml to validate
    :param type: kind of yaml
    :param local: files are local
    :return:
    """
    assert type in ["alignment", "hmmcopy", "variant_calling",
                    "breakpoint_calling", "count_haps", "merge_cell_bams",
                    "split_wgs_bams", "annotation"]
    yaml = load(yaml)
    yaml = make_parsed_format_consistent(yaml)

    input_files = list(filter(could_be_file , list(flatten(yaml).values())))

    output_files = []

    #count haps,
    if type == "count_haps":
        output_files = [yaml["haplotypes"]]

    assert check_local_input_files(input_files)
    assert check_local_output_file(input_files)

    #special check for alignment types
    if type == "alignment":
        alignment_type_example = load("alignment_yaml_types.yaml")

            assert validate_alignment_types(yaml[cell], alignment_type_example["CELL"])

if __name__=="__main__":
    validate_yaml(sys.argv[1], sys.argv[2])
