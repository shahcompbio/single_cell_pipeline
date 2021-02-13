import yaml
from single_cell.utils.validator import validate


def load_split_wgs_input(input_yaml):
    yamldata = load_yaml(input_yaml)

    validate.validate_split_wgs_bam(yamldata)

    wgs_bams = yamldata['normal']

    return wgs_bams['bam']


def load_merge_cell_bams(input_yaml):
    yamldata = load_yaml(input_yaml)

    validate.validate_merge_cell_bams(yamldata)

    cell_bams = yamldata['cell_bams']

    cell_bams = {cell_id: cell_bams[cell_id]['bam'] for cell_id in cell_bams}

    return cell_bams


def load_infer_haps_input(input_yaml):
    yamldata = load_yaml(input_yaml)

    validate.validate_infer_haps(yamldata)

    normal = yamldata['normal']

    if 'bam' in normal:
        normal = normal['bam']
    else:
        normal = {v: normal[v]['bam'] for v in normal}

    return normal


def load_count_haps_input(input_yaml):
    yamldata = load_yaml(input_yaml)

    validate.validate_count_haps(yamldata)

    haplotypes = yamldata['haplotypes']

    tumours = yamldata['tumour']

    tumours = {v: tumours[v]['bam'] for v in tumours}

    return haplotypes, tumours


def load_breakpoint_calling_input(input_yaml):
    yamldata = load_yaml(input_yaml)

    validate.validate_breakpoint_calling(yamldata)

    normal = yamldata['normal']

    if 'bam' in normal:
        normal = normal['bam']
    else:
        normal = {v: normal[v]['bam'] for v in normal}

    tumours = yamldata['tumour']

    tumours = {v: tumours[v]['bam'] for v in tumours}

    return normal, tumours


def load_config(args):
    return load_yaml(args["config_file"])


def load_variant_calling_input(input_yaml):
    yamldata = load_yaml(input_yaml)

    validate.validate_variant_calling(yamldata)

    normals = yamldata['normal']

    tumours = yamldata['tumour']

    normals = {v: normals[v]['bam'] for v in normals}
    tumours = {v: tumours[v]['bam'] for v in tumours}

    assert normals.keys() == tumours.keys()

    return normals, tumours


def load_germline_data(input_yaml):
    yamldata = load_yaml(input_yaml)

    validate.validate_germline_calling(yamldata)

    normals = yamldata['normal']
    normals = {v: normals[v]['bam'] for v in normals}

    return normals


def load_variant_counting_input(input_yaml):
    yamldata = load_yaml(input_yaml)

    validate.validate_snv_genotyping(yamldata)

    vcf_files = yamldata['vcf_files']

    strelka_vcf_data = {}
    museq_vcf_data = {}
    for sample, sampledata in vcf_files.items():
        for library, library_data in sampledata.items():
            strelka_vcf_data[(sample, library)] = library_data['strelka_snv_vcf']
            museq_vcf_data[(sample, library)] = library_data['museq_vcf']

    cells_data = yamldata['tumour_cells']

    sample_library = []
    cells_data_out = {}
    for sample, sampledata in cells_data.items():
        for library, library_data in sampledata.items():
            sample_library.append({'sample_id': sample, 'library_id': library})
            for cell, cell_data in library_data.items():
                cells_data_out[(sample, library, cell)] = cell_data['bam']

    return strelka_vcf_data, museq_vcf_data, cells_data_out, sample_library


def load_sv_genotyper_input(input_yaml):
    yamldata = load_yaml(input_yaml)

    validate.validate_sv_genotyping(yamldata)

    sv_calls = yamldata['sv_calls']
    lumpy_csvs = {}
    destruct_csvs = {}
    for sample, sampledata in sv_calls.items():
        for library, library_data in sampledata.items():
            lumpy_csvs[(sample, library)] = library_data['lumpy']
            destruct_csvs[(sample, library)] = library_data['destruct']

    cells_data = yamldata['tumour_cells']

    cells_data_out = {}
    for sample, sampledata in cells_data.items():
        for library, library_data in sampledata.items():
            for cell, cell_data in library_data.items():
                cells_data_out[(sample, library, cell)] = cell_data['bam']

    return lumpy_csvs, destruct_csvs, cells_data_out


def load_yaml(path):
    try:
        with open(path) as infile:
            data = yaml.safe_load(infile)

    except IOError:
        raise Exception(
            'Unable to open file: {0}'.format(path))
    return data


def load_qc_input(path):
    data = {}
    yaml = load_yaml(path)
    for patient, patient_data in yaml.items():
        for sample, sample_data in patient_data.items():
            for library, library_data in sample_data.items():
                data[(patient, sample, library)] = {
                    data_label: data for data_label, data in library_data.items()
                }
    return data


def get_sample_info(fastqs_file):
    """
    load yaml and remove some extra info to reduce size
    """

    data = load_yaml(fastqs_file)

    validate.validate_sample_info(data)

    cells = data.keys()

    for cell in cells:
        data[cell]["cell_call"] = data[cell]["pick_met"]
        data[cell]["experimental_condition"] = data[cell]["condition"]
        if "fastqs" in data[cell]:
            del data[cell]["fastqs"]
        if "bam" in data[cell]:
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
        assert "bam" in data[cell], \
            "couldnt extract bam file paths from yaml input for cell: {}".format(cell)

    bam_filenames = {cell: data[cell]["bam"] for cell in data.keys()}

    return bam_filenames


def get_fastqs(fastqs_file):
    data = load_yaml(fastqs_file)

    validate.validate_alignment_fastqs(data)

    for cell in data.keys():
        assert "fastqs" in data[cell], \
            "couldnt extract fastq file paths from yaml input for cell: {}".format(cell)

    fastq_1_filenames = dict()
    fastq_2_filenames = dict()
    for cell in data.keys():
        fastqs = data[cell]["fastqs"]

        for lane, paths in fastqs.items():
            fastq_1_filenames[(cell, lane)] = paths["fastq_1"]
            fastq_2_filenames[(cell, lane)] = paths["fastq_2"]

    return fastq_1_filenames, fastq_2_filenames


def load_cohort_qc_inputs(yaml_file):
    yaml_data = load_yaml(yaml_file)

    mafs = {}
    hmmcopy_filenames = {}
    hmmcopy_samples = {}

    for sample, data in yaml_data["hmmcopy"].items():
        for library, data in data.items():
            hmmcopy_filenames[(sample, library)] = data["hmmcopy_reads"]
            hmmcopy_samples[(sample, library)] = data["sample_id"]

    for sample, data in yaml_data["mafs"].items():
        mafs[(sample)] = data

    return mafs, hmmcopy_filenames, hmmcopy_samples

