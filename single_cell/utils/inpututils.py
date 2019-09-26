import yaml


def load_config(args):
    return load_yaml(args["config_file"])


def load_yaml_section(data, section_name):
    if data.get(section_name):
        assert len(data[section_name]) == 1
        section_id = data[section_name].keys()[0]
        section_data = data[section_name][section_id]
        if 'bam' in section_data:
            section_data = section_data['bam']
        else:
            section_data = {cell_id: bamdata['bam'] for cell_id, bamdata in section_data.items()}
    else:
        section_id = None
        section_data = None
    return section_id, section_data


def load_pseudowgs_input(inputs_file):
    data = load_yaml(inputs_file)

    normal_wgs_id, normal_wgs = load_yaml_section(data, 'normal_wgs')
    tumour_wgs_id, tumour_wgs = load_yaml_section(data, 'tumour_wgs')

    tumour_cells_id, tumour_cells = load_yaml_section(data, 'tumour_cells')
    normal_cells_id, normal_cells = load_yaml_section(data, 'normal_cells')

    parsed_data = dict(
        tumour_wgs_id=tumour_wgs_id, tumour_wgs=tumour_wgs,
        normal_wgs_id=normal_wgs_id, normal_wgs=normal_wgs,
        tumour_cells_id=tumour_cells_id, tumour_cells=tumour_cells,
        normal_cells_id=normal_cells_id, normal_cells=normal_cells)

    return parsed_data


def load_yaml(path):
    try:
        with open(path) as infile:
            data = yaml.safe_load(infile)

    except IOError:
        raise Exception(
            'Unable to open file: {0}'.format(path))
    return data


def get_fastqs(fastqs_file):
    data = load_yaml(fastqs_file)

    for cell in data.keys():
        assert "fastqs" in data[
            cell], "couldnt extract fastq file paths from yaml input for cell: {}".format(cell)

    fastq_1_filenames = dict()
    fastq_2_filenames = dict()
    for cell in data.keys():
        fastqs = data[cell]["fastqs"]

        for lane, paths in fastqs.items():
            fastq_1_filenames[(cell, lane)] = paths["fastq_1"]
            fastq_2_filenames[(cell, lane)] = paths["fastq_2"]

    return fastq_1_filenames, fastq_2_filenames


def get_trim_info(fastqs_file):
    data = load_yaml(fastqs_file)

    for cell in data.keys():
        assert "fastqs" in data[
            cell], "couldnt extract fastq file paths from yaml input for cell: {}".format(cell)

    seqinfo = dict()
    for cell in data.keys():
        fastqs = data[cell]["fastqs"]

        for lane, paths in fastqs.items():
            if 'trim' in paths:
                seqinfo[(cell, lane)] = paths["trim"]
            elif 'sequencing_instrument' in paths:
                DeprecationWarning("sequencing instrument value is deprecated "
                                   "and will be removed with v0.2.8")
                if paths["sequencing_instrument"] == "N550":
                    trim = False
                else:
                    trim = True
                seqinfo[(cell, lane)] = trim
            else:
                raise Exception(
                    "trim flag missing in cell: {}".format(cell))

    return seqinfo


def get_center_info(fastqs_file):
    data = load_yaml(fastqs_file)

    for cell in data.keys():
        assert "fastqs" in data[
            cell], "couldnt extract fastq file paths from yaml input for cell: {}".format(cell)

    seqinfo = dict()
    for cell in data.keys():
        fastqs = data[cell]["fastqs"]

        for lane, paths in fastqs.items():

            if "sequencing_center" not in paths:
                raise Exception(
                    "sequencing_center key missing in cell: {}".format(cell))
            seqinfo[(cell, lane)] = paths["sequencing_center"]

    return seqinfo


def get_sample_info(fastqs_file):
    """
    load yaml and remove some extra info to reduce size
    """

    data = load_yaml(fastqs_file)

    cells = data.keys()

    for cell in cells:
        data[cell]["cell_call"] = data[cell]["pick_met"]
        data[cell]["experimental_condition"] = data[cell]["condition"]
        if "fastqs" in data[cell]:
            del data[cell]["fastqs"]
        del data[cell]["bam"]
        del data[cell]["pick_met"]
        del data[cell]["condition"]

    return data


def get_samples(fastqs_file):
    data = load_yaml(fastqs_file)

    return list(data.keys())


def get_bams(fastqs_file):
    data = load_yaml(fastqs_file)

    for cell in data.keys():
        assert "bam" in data[
            cell], "couldnt extract bam file paths from yaml input for cell: {}".format(cell)

    bam_filenames = {cell: data[cell]["bam"] for cell in data.keys()}
    bai_filenames = {cell: data[cell]["bam"] + ".bai" for cell in data.keys()}

    return bam_filenames, bai_filenames


def load_cell_data(yamldata, key):
    for sample_id, sample_data in yamldata[key].items():
        for library_id, library_data in sample_data.items():
            for cell_id, cell_data in library_data.items():
                yield sample_id, library_id, cell_id, cell_data['bam']


def load_tumour_data(yamldata):
    sample_data = {}

    for sample_id, library_id, cell_id, cell_bam in load_cell_data(yamldata, 'tumour_cells'):
        sample_data[(sample_id, library_id, cell_id)] = cell_bam

    return sample_data


def load_normal_data(yamldata):
    libraries = set()
    if 'normal_wgs' in yamldata:
        assert len(yamldata['normal_wgs'].keys()) == 1
        sample_id = list(yamldata['normal_wgs'].keys())[0]
        assert len(yamldata['normal_wgs'][sample_id].keys()) == 1
        library_id = list(yamldata['normal_wgs'][sample_id].keys())[0]
        libraries.add(library_id)
        cell_bams = yamldata['normal_wgs'][sample_id][library_id]['bam']
    else:
        cell_bams = {}

        if not len(yamldata['normal_cells'].keys()) == 1:
            raise Exception("Pipeline does not support multiple normal samples")

        for sample_id, library_id, cell_id, cell_bam in load_cell_data(yamldata, 'normal_cells'):
            libraries.add(library_id)
            if cell_id in cell_bams:
                raise Exception("non unique cell id {} encountered".format(cell_id))
            cell_bams[cell_id] = cell_bam

    return sample_id, sorted(libraries), cell_bams
