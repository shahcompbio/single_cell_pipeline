import os
import pypeliner
import pypeliner.managed as mgd

from single_cell.utils import helpers
from single_cell.utils import refgenome
from workflows import split_bams
from workflows import merge_bams
import infer_haps
import variant_calling
import single_cell.workflows.split_bams.tasks


def get_bams(inputs_file):
    data = helpers.load_yaml(inputs_file)

    tumour_bams = {}
    for sample_id in data['tumour'].keys():
        for cell in data['tumour'][sample_id].keys():
            if 'bam' not in data['tumour'][sample_id][cell]:
                raise Exception('couldnt extract bam file paths from yaml input for cell: {}'.format(cell))
            tumour_bams[(sample_id, cell)] = data['tumour'][sample_id][cell]['bam']

    if 'bam' in data['normal']:
        normal_bams = data['normal']['bam']
    else:
        normal_bams = {}
        for cell in data['tumour'].keys():
            if 'bam' not in data['tumour'][cell]:
                raise Exception('couldnt extract bam file paths from yaml input for cell: {}'.format(cell))
            normal_bams[cell] = data['tumour'][cell]['bam']

    return tumour_bams, normal_bams


def multi_sample_workflow(args):
    tumour_cell_bams, normal_bam = get_bams(args['input_yaml'])

    workflow = create_multi_sample_workflow(
        normal_bam,
        tumour_cell_bams,
        args['out_dir'],
        helpers.load_config(args),
    )

    return workflow


def create_multi_sample_workflow(
        normal_wgs_bam,
        tumour_cell_bams,
        results_dir,
        config,
):
    """ Multiple sample pseudobulk workflow. """

    baseimage = config['multi_sample']['docker']['single_cell_pipeline']
    ctx = {'mem_retry_increment': 2, 'ncpus': 1}

    raw_data_dir = os.path.join(results_dir, 'raw')

    normal_region_bam_template = os.path.join(raw_data_dir, 'normal_{region}.bam')
    tumour_region_bam_template = os.path.join(raw_data_dir, '{sample_id}_{region}.bam')
    normal_seqdata_file = os.path.join(raw_data_dir, 'normal_seqdata.h5')
    tumour_cell_seqdata_template = os.path.join(raw_data_dir, '{sample_id}_{cell_id}_seqdata.h5')

    variant_calling_raw_data_template = os.path.join(raw_data_dir, '{sample_id}_variant_calling')
    destruct_raw_data_template = os.path.join(raw_data_dir, '{sample_id}_destruct')

    museq_vcf_template = os.path.join(results_dir, '{sample_id}_museq.vcf.gz')
    strelka_snv_template = os.path.join(results_dir, '{sample_id}_strelka_snv.vcf.gz')
    strelka_indel_template = os.path.join(results_dir, '{sample_id}_strelka_indel.vcf.gz')
    snv_annotations_template = os.path.join(results_dir, '{sample_id}_snv_annotations.h5')
    snv_counts_template = os.path.join(results_dir, '{sample_id}_snv_counts.h5')
    haplotypes_file = os.path.join(results_dir, 'haplotypes.tsv')
    allele_counts_template = os.path.join(results_dir, '{sample_id}_allele_counts.csv')
    breakpoints_template = os.path.join(results_dir, '{sample_id}_destruct.h5')

    snv_calling_info_template = os.path.join(results_dir, '{sample_id}_snv_calling_info.yaml')
    snv_counting_info_template = os.path.join(results_dir, '{sample_id}_snv_counting_info.yaml')
    multisample_info_filename = os.path.join(results_dir, 'multisample_info.yaml')

    regions = refgenome.get_split_regions(config["split_bam"]["split_size"])

    workflow = pypeliner.workflow.Workflow(default_ctx=ctx)

    workflow.set_filenames('normal_regions.bam', 'region', template=normal_region_bam_template)
    workflow.set_filenames('tumour_cells.bam', 'sample_id', 'cell_id', fnames=tumour_cell_bams)
    workflow.set_filenames('tumour_regions.bam', 'sample_id', 'region', template=tumour_region_bam_template)

    workflow.set_filenames('museq.vcf', 'sample_id', template=museq_vcf_template)
    workflow.set_filenames('strelka_snv.vcf', 'sample_id', template=strelka_snv_template)
    workflow.set_filenames('strelka_indel.vcf', 'sample_id', template=strelka_indel_template)
    workflow.set_filenames('snv_annotations.h5', 'sample_id', template=snv_annotations_template)
    workflow.set_filenames('snv_counts.h5', 'sample_id', template=snv_counts_template)
    workflow.set_filenames('tumour_cell_seqdata.h5', 'sample_id', 'cell_id', template=tumour_cell_seqdata_template)
    workflow.set_filenames('allele_counts.csv', 'sample_id', template=allele_counts_template)
    workflow.set_filenames('breakpoints.h5', 'sample_id', template=breakpoints_template)

    workflow.set_filenames('snv_calling_info.yaml', 'sample_id', template=snv_calling_info_template)
    workflow.set_filenames('snv_counting_info.yaml', 'sample_id', template=snv_counting_info_template)

    workflow.setobj(
        obj=mgd.OutputChunks('region'),
        value=regions,
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'cell_id'),
        value=tumour_cell_bams.keys(),
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'region'),
        axes=('sample_id',),
        value=regions,
    )

    if isinstance(normal_wgs_bam, dict):
        workflow.set_filenames('normal_cells.bam', 'cell_id', fnames=normal_wgs_bam)
        workflow.subworkflow(
            name="split_normal",
            func=merge_bams.create_merge_bams_workflow,
            args=(
                mgd.InputFile('normal_cells.bam', 'cell_id', extensions=['.bai']),
                mgd.OutputFile('normal_regions.bam', 'region', axes_origin=[], extensions=['.bai']),
                regions,
                config['merge_bams'],
            )
        )
    else:
        workflow.subworkflow(
            name="split_normal",
            func=split_bams.create_split_workflow,
            args=(
                mgd.InputFile(normal_wgs_bam, extensions=['.bai']),
                mgd.OutputFile('normal_regions.bam', 'region', extensions=['.bai'], axes_origin=[]),
                pypeliner.managed.TempInputObj('region'),
                config['split_bam'],
            ),
            kwargs={"by_reads": False}
        )

    workflow.subworkflow(
        name="split_merge_tumour",
        axes=('sample_id',),
        func=merge_bams.create_merge_bams_workflow,
        args=(
            mgd.InputFile('tumour_cells.bam', 'sample_id', 'cell_id', extensions=['.bai']),
            mgd.OutputFile('tumour_regions.bam', 'sample_id', 'region', axes_origin=[], extensions=['.bai']),
            regions,
            config['merge_bams'],
        )
    )

    workflow.subworkflow(
        name='variant_calling',
        func=variant_calling.create_variant_calling_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile('tumour_cells.bam', 'sample_id', 'cell_id', extensions=['.bai']),
            mgd.InputFile('tumour_regions.bam', 'sample_id', 'region', extensions=['.bai']),
            mgd.InputFile('normal_regions.bam', 'region', extensions=['.bai']),
            mgd.OutputFile('museq.vcf', 'sample_id'),
            mgd.OutputFile('strelka_snv.vcf', 'sample_id'),
            mgd.OutputFile('strelka_indel.vcf', 'sample_id'),
            mgd.OutputFile('snv_annotations.h5', 'sample_id'),
            mgd.OutputFile('snv_calling_info.yaml', 'sample_id'),
            config['variant_calling'],
            mgd.Template(variant_calling_raw_data_template, 'sample_id'),
        ),
    )

    workflow.transform(
        name='merge_museq_snvs',
        func='biowrappers.components.io.vcf.tasks.concatenate_vcf',
        args=(
            mgd.InputFile('museq.vcf', 'sample_id', axes_origin=[]),
            mgd.TempOutputFile('museq.vcf.gz', extensions=['.tbi', '.csi']),
        ),
    )

    workflow.transform(
        name='merge_strelka_snvs',
        func='biowrappers.components.io.vcf.tasks.concatenate_vcf',
        args=(
            mgd.InputFile('strelka_snv.vcf', 'sample_id', axes_origin=[]),
            mgd.TempOutputFile('strelka_snv.vcf.gz', extensions=['.tbi', '.csi']),
        ),
    )

    workflow.subworkflow(
        name='variant_counting',
        func=variant_calling.create_variant_counting_workflow,
        axes=('sample_id',),
        args=(
            [
                mgd.TempInputFile('museq.vcf.gz', extensions=['.tbi', '.csi']),
                mgd.TempInputFile('strelka_snv.vcf.gz', extensions=['.tbi', '.csi']),
            ],
            mgd.InputFile('tumour_cells.bam', 'sample_id', 'cell_id', extensions=['.bai']),
            mgd.OutputFile('snv_counts.h5', 'sample_id'),
            mgd.OutputFile('snv_counting_info.yaml', 'sample_id'),
            config['variant_calling'],
        ),
    )

    workflow.subworkflow(
        name='infer_haps_from_bulk_normal',
        func=infer_haps.infer_haps_from_bulk_normal,
        args=(
            mgd.InputFile(normal_wgs_bam, extensions=['.bai']),
            mgd.OutputFile(normal_seqdata_file),
            mgd.OutputFile(haplotypes_file),
            config['infer_haps'],
        ),
    )

    workflow.subworkflow(
        name='extract_allele_readcounts',
        func=infer_haps.extract_allele_readcounts,
        axes=('sample_id',),
        args=(
            mgd.InputFile(haplotypes_file),
            mgd.InputFile('tumour_cells.bam', 'sample_id', 'cell_id', extensions=['.bai']),
            mgd.OutputFile('tumour_cell_seqdata.h5', 'sample_id', 'cell_id', axes_origin=[]),
            mgd.OutputFile('allele_counts.csv', 'sample_id', template=allele_counts_template),
            config['infer_haps'],
        ),
    )

    destruct_config = config['breakpoint_calling'].get('destruct_config', {})
    destruct_ref_data_dir = config['breakpoint_calling']['ref_data_directory']

    workflow.subworkflow(
        name='destruct',
        func='biowrappers.components.breakpoint_calling.destruct.destruct_pipeline',
        axes=('sample_id',),
        args=(
            mgd.InputFile(normal_wgs_bam),
            mgd.InputFile('tumour_cells.bam', 'sample_id', 'cell_id', extensions=['.bai']),
            destruct_config,
            destruct_ref_data_dir,
            mgd.OutputFile('breakpoints.h5', 'sample_id'),
            mgd.Template(destruct_raw_data_template, 'sample_id'),
        ),
    )

    workflow.transform(
        name="get_allele_count_filename",
        func='single_cell.utils.helpers.resolve_template',
        ctx=dict(docker_image=baseimage, **ctx),
        ret=pypeliner.managed.TempOutputObj('allele_counts'),
        args=(
            pypeliner.managed.TempInputObj('sample_id'),
            allele_counts_template,
            'sample_id'
        )
    )

    workflow.transform(
        name="get_breakpoint_filename",
        func='single_cell.utils.helpers.resolve_template',
        ctx=dict(docker_image=baseimage, **ctx),
        ret=pypeliner.managed.TempOutputObj('breakpoints'),
        args=(
            pypeliner.managed.TempInputObj('sample_id'),
            breakpoints_template,
            'sample_id'
        )
    )

    if isinstance(normal_wgs_bam, dict):
        normal_wgs_bam = {
            cell: helpers.format_file_yaml(bam) for cell, bam in normal_wgs_bam.iteritems()
        }
    else:
        normal_wgs_bam = helpers.format_file_yaml(normal_wgs_bam)

    tumour_cell_bams = {
        ','.join(k): helpers.format_file_yaml(bam) for k, bam in tumour_cell_bams.iteritems()
    }

    metadata = {
                   'name': 'multi_sample_pseudobulk',
                   'version': single_cell.__version__,
                   'output_datasets': None,
                   'input_datasets': {
                       'normal_bam': normal_wgs_bam,
                       'tumour_cell_bams': tumour_cell_bams,
                   },
                   'results': {
                       'haplotypes': helpers.format_file_yaml(haplotypes_file),
                       'allele_counts': mgd.TempInputObj('allele_counts'),
                       'breakpoints': mgd.TempInputObj('breakpoints'),
                   },
               },

    workflow.transform(
        name='generate_meta_yaml',
        func="single_cell.utils.helpers.write_to_yaml",
        ctx=dict(docker_image=baseimage, **ctx),
        args=(
            mgd.OutputFile(multisample_info_filename),
            metadata
        )
    )

    return workflow
