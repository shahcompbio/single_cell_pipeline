from single_cell.utils.validator import utils


def validate_alignment_fastqs(data):
    for sample, sample_data in data.items():
        for lane, lane_data in sample_data['fastqs'].items():
            if not utils.get(lane_data, 'fastq_1') or not utils.get(lane_data, 'fastq_2'):
                raise utils.MissingInput()
            utils.check_data_type(['sequencing_center', 'sequencing_instrument'], str, lane_data)

            # this is not the platform. that is hardcoded in pipeline.
            # utils.check_sequencing_instrument_type(utils.get(lane_data, 'sequencing_instrument'))


def validate_sample_info(yamldata):
    for cell in yamldata:
        celldata = yamldata[cell]

        utils.check_data_type(['column', 'img_col', 'row'], int, celldata)
        utils.check_data_type(['condition', 'pick_met', 'index_i5', 'index_i7'], str, celldata)

        utils.check_barcodes(utils.get(celldata, 'primer_i5'))
        utils.check_barcodes(utils.get(celldata, 'primer_i7'))

        if not utils.get(celldata, 'index_i5').startswith('i5-'):
            raise utils.DLPIndexError()
        if not utils.get(celldata, 'index_i7').startswith('i7-'):
            raise utils.DLPIndexError()


def validate_hmmcopy_bams(yamldata):
    for cell, celldata in yamldata.items():
        utils.check_data_type(['bam'], str, celldata)


def validate_annotation(yamldata):
    utils.check_data_type(
        ['hmmcopy_metrics', 'hmmcopy_reads', 'alignment_metrics', 'gc_metrics', 'segs_pdf_tar'],
        str,
        yamldata
    )


def validate_merge_cell_bams(yamldata):
    utils.check_cells_data(utils.get(yamldata, 'cell_bams'))


def validate_split_wgs_bam(yamldata):
    data = utils.get(yamldata, 'normal')
    utils.check_data_type(['bam'], str, data)


def validate_variant_calling(yamldata):
    normals = yamldata['normal']
    for region in normals:
        utils.check_data_type(['bam'], str, normals[region])
        utils.check_genomic_regions(region)

    tumours = yamldata['tumour']
    for region in tumours:
        utils.check_data_type(['bam'], str, tumours[region])
        utils.check_genomic_regions(region)


def validate_germline_calling(yamldata):
    utils.check_normal_data(utils.get(yamldata, 'normal'))


def validate_infer_haps(yamldata):
    utils.check_normal_data(utils.get(yamldata, 'normal'))


def validate_count_haps(yamldata):
    utils.check_cells_data(utils.get(yamldata, 'tumour'))
    utils.check_data_type(['haplotypes'], str, yamldata)


def validate_breakpoint_calling(yamldata):
    utils.check_normal_data(utils.get(yamldata, 'normal'))
    utils.check_cells_data(utils.get(yamldata, 'tumour'))


def validate_snv_genotyping(yamldata):
    tumour_cells = utils.get(yamldata, 'tumour_cells')
    for sample in tumour_cells:
        for library in tumour_cells[sample]:
            utils.check_cells_data(tumour_cells[sample][library])

    vcf_files = utils.get(yamldata, 'vcf_files')
    for sample in vcf_files:
        for library in vcf_files[sample]:
            utils.check_data_type(['museq_vcf', 'strelka_snv_vcf'], str, vcf_files[sample][library])


def validate_sv_genotyping(yamldata):
    pass
