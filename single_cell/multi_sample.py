import os
import pypeliner
import pypeliner.managed as mgd

from single_cell.utils import helpers
from single_cell.utils import refgenome
import infer_haps
import variant_calling
import single_cell.workflows.split_bams.tasks


def get_bams(inputs_file):
    data = helpers.load_yaml(inputs_file)

    bam_filenames = {}
    for sample_id in data['tumour'].keys():
        for cell in data['tumour'][sample_id].keys():
            if 'bam' not in data['tumour'][sample_id][cell]:
                raise Exception('couldnt extract bam file paths from yaml input for cell: {}'.format(cell))
            bam_filenames[(sample_id, cell)] = data['tumour'][sample_id][cell]['bam']

    normal_bam = data['normal']['bam']

    return bam_filenames, normal_bam


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

    ctx = {'mem_retry_increment': 2, 'ncpus': 1}
    docker_ctx = helpers.get_container_ctx(config['containers'], 'single_cell_pipeline')
    ctx.update(docker_ctx)

    raw_data_dir = os.path.join(results_dir, 'raw')

    normal_region_bam_template = os.path.join(raw_data_dir, 'normal_{region}.bam')
    tumour_region_bam_template = os.path.join(raw_data_dir, '{sample_id}_{region}.bam')
    normal_seqdata_file = os.path.join(raw_data_dir, 'normal_seqdata.h5')
    tumour_cell_seqdata_template = os.path.join(raw_data_dir, '{sample_id}_{cell_id}_seqdata.h5')

    museq_vcf_template = os.path.join(results_dir, '{sample_id}_museq.vcf.gz')
    strelka_snv_template = os.path.join(results_dir, '{sample_id}_strelka_snv.vcf.gz')
    strelka_indel_template = os.path.join(results_dir, '{sample_id}_strelka_indel.vcf.gz')
    snv_annotations_template = os.path.join(results_dir, '{sample_id}_snv_annotations.h5')
    snv_meta_template = os.path.join(results_dir, '{sample_id}_snv_meta.yaml')
    snv_counts_template = os.path.join(results_dir, '{sample_id}_snv_counts.h5')
    haplotypes_file = os.path.join(results_dir, 'haplotypes.tsv')
    allele_counts_template = os.path.join(results_dir, '{sample_id}_allele_counts.tsv')
    breakpoints_template = os.path.join(results_dir, '{sample_id}.h5')

    regions = refgenome.get_split_regions(config["split_size"])

    workflow = pypeliner.workflow.Workflow(default_ctx=ctx)

    workflow.set_filenames('normal_regions.bam', 'region', template=normal_region_bam_template)
    workflow.set_filenames('tumour_cells.bam', 'sample_id', 'cell_id', fnames=tumour_cell_bams)
    workflow.set_filenames('tumour_regions.bam', 'sample_id', 'region', template=tumour_region_bam_template)

    workflow.set_filenames('museq.vcf', 'sample_id', template=museq_vcf_template)
    workflow.set_filenames('strelka_snv.vcf', 'sample_id', template=strelka_snv_template)
    workflow.set_filenames('strelka_indel.vcf', 'sample_id', template=strelka_indel_template)
    workflow.set_filenames('snv_annotations.h5', 'sample_id', template=snv_annotations_template)
    workflow.set_filenames('snv_meta.yaml', 'sample_id', template=snv_meta_template)
    workflow.set_filenames('snv_counts.h5', 'sample_id', template=snv_counts_template)
    workflow.set_filenames('tumour_cell_seqdata.h5', 'sample_id', 'cell_id', template=tumour_cell_seqdata_template)
    workflow.set_filenames('breakpoints.h5', 'sample_id', template=breakpoints_template)

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

    workflow.transform(
        name='split_normal',
        func='single_cell.workflows.split_bams.tasks.split_bam_file_one_job',
        args=(
            mgd.InputFile(normal_wgs_bam, extensions=['.bai']),
            mgd.OutputFile('normal_regions.bam', 'region', extensions=['.bai'], axes_origin=[]),
            mgd.InputChunks('region'),
            helpers.get_container_ctx(config['containers'], 'samtools'),
        ),
        kwargs={
            'ncores': config['max_cores'],
        },
    )

    workflow.subworkflow(
        name='split_merge_tumour',
        func='single_cell.workflows.merge_bams.create_cell_region_merge_workflow',
        axes=('sample_id',),
        args=(
            mgd.InputFile('tumour_cells.bam', 'sample_id', 'cell_id', extensions=['.bai']),
            mgd.OutputFile('tumour_regions.bam', 'sample_id', 'region', axes_origin=[], extensions=['.bai']),
            regions,
            config,
        ),
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
            mgd.OutputFile('snv_meta.yaml', 'sample_id'),
            config,
            raw_data_dir,
        ),
    )

    workflow.transform(
        name='merge_museq_snvs',
        func='biowrappers.components.io.vcf.tasks.concatenate_vcf',
        args=(
            mgd.InputFile('museq.vcf', 'sample_id', axes_origin=[]),
            mgd.TempOutputFile('museq.vcf'),
        ),
    )

    workflow.transform(
        name='merge_strelka_snvs',
        func='biowrappers.components.io.vcf.tasks.concatenate_vcf',
        args=(
            mgd.InputFile('strelka_snv.vcf', 'sample_id', axes_origin=[]),
            mgd.TempOutputFile('strelka_snv.vcf'),
        ),
    )

    workflow.subworkflow(
        name='variant_counting',
        func=variant_calling.create_variant_counting_workflow,
        axes=('sample_id',),
        args=(
            [
                mgd.TempInputFile('museq.vcf'),
                mgd.TempInputFile('strelka_snv.vcf'),
            ],
            mgd.InputFile('tumour_cells.bam', 'sample_id', 'cell_id', extensions=['.bai']),
            mgd.OutputFile('snv_counts.h5', 'sample_id'),
            config,
        ),
    )

    workflow.subworkflow(
        name='infer_haps_from_bulk_normal',
        func=infer_haps.infer_haps_from_bulk_normal,
        args=(
            mgd.InputFile(normal_wgs_bam, extensions=['.bai']),
            mgd.OutputFile(normal_seqdata_file),
            mgd.OutputFile(haplotypes_file),
            config,
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
            mgd.OutputFile('allele_counts.h5', 'sample_id', template=allele_counts_template),
            config,
        ),
    )

    destruct_config = config.get('destruct_config', {})
    destruct_ref_data_dir = config['destruct_ref_data_dir']

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
            raw_data_dir,
        ),
    )

    return workflow

