'''
Created on Feb 22, 2018

@author: dgrewal
'''

import os

import pypeliner.managed as mgd
from single_cell.utils import helpers
from single_cell.utils import refgenome
from single_cell.workflows import merge_bams
from single_cell.workflows import mutationseq
from single_cell.workflows import split_bams
from single_cell.workflows import strelka

import pypeliner

default_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']


def create_snv_allele_counts_for_vcf_targets_workflow(
        bam_files,
        vcf_file,
        out_file,
        memory_cfg,
        count_duplicates=False,
        min_bqual=0,
        min_mqual=0,
        table_name='snv_allele_counts',
        vcf_to_bam_chrom_map=None,
):
    ctx = {
        'mem': memory_cfg['low'], 'num_retry': 3, 'mem_retry_increment': 2, 'ncpus': 1,
        'disk_retry_increment': 50,
    }
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=list(bam_files.keys()),
    )

    workflow.transform(
        name='get_snv_allele_counts_for_vcf_targets',
        axes=('cell_id',),
        func="biowrappers.components.variant_calling.snv_allele_counts.tasks.get_snv_allele_counts_for_vcf_targets",
        args=(
            mgd.InputFile('tumour.bam', 'cell_id', fnames=bam_files, extensions=['.bai']),
            mgd.InputFile(vcf_file),
            mgd.TempOutputFile('counts.h5', 'cell_id'),
            table_name,
        ),
        kwargs={
            'count_duplicates': count_duplicates,
            'min_bqual': min_bqual,
            'min_mqual': min_mqual,
            'vcf_to_bam_chrom_map': vcf_to_bam_chrom_map,
            'cell_id': mgd.Instance('cell_id'),
            'report_zero_count_positions': False,
        }
    )

    workflow.transform(
        name='merge_snv_allele_counts',
        ctx={'mem': memory_cfg['high'], 'disk': 20},
        func="biowrappers.components.io.hdf5.tasks.concatenate_tables",
        args=(
            mgd.TempInputFile('counts.h5', 'cell_id'),
            mgd.OutputFile(out_file),
        ),
        kwargs={
            'in_memory': False,
        },
    )

    return workflow


def museq_callback(record):
    return record.INFO['PR']


def strelka_snv_callback(record):
    return record.INFO['QSS']


def variant_calling_workflow(args):
    config = helpers.load_config(args)
    config = config['variant_calling']

    meta_yaml = os.path.join(args['out_dir'], 'info.yaml')

    data = helpers.load_pseudowgs_input(args['input_yaml'])
    tumour_bams = data['tumour_wgs']
    normal_bams = data['normal_wgs']
    tumour_cells = data['tumour_cells']

    varcalls_dir = os.path.join(args['out_dir'], 'results',
                                'variant_calling')

    museq_vcf = os.path.join(varcalls_dir, 'museq_snv.vcf.gz')
    strelka_snv_vcf = os.path.join(varcalls_dir, 'strelka_snv.vcf.gz')
    strelka_indel_vcf = os.path.join(varcalls_dir, 'strelka_indel.vcf.gz')
    snv_h5 = os.path.join(varcalls_dir, 'snv_annotations.h5')
    raw_data_dir = os.path.join(varcalls_dir, 'raw')

    baseimage = config['docker']['single_cell_pipeline']

    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1,
           'mem': config["memory"]['low'], 'docker_image': baseimage}
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    if isinstance(normal_bams, dict) and isinstance(tumour_bams, dict):
        assert list(normal_bams.keys()) == list(
            tumour_bams.keys()), 'keys for tumour and normal bams should be the same'
        workflow.setobj(
            obj=mgd.OutputChunks('region'),
            value=list(normal_bams.keys()),
        )
        workflow.set_filenames('normal_split.bam', 'normal_split', fnames=normal_bams)
        workflow.set_filenames('tumour_split.bam', 'normal_split', fnames=tumour_bams)
    else:
        workflow.transform(
            name="get_regions",
            func="single_cell.utils.pysamutils.get_regions_from_reference",
            ret=pypeliner.managed.OutputChunks('region'),
            args=(
                config["ref_genome"],
                config["split_size"],
                config["chromosomes"],
            )
        )
        assert '{region}' in normal_bams, 'only supports a list of files or a template on regions'
        workflow.set_filenames('normal_split.bam', 'region', template=normal_bams)
        assert '{region}' in tumour_bams, 'only supports a list of files or a template on regions'
        workflow.set_filenames('tumour_split.bam', 'region', template=normal_bams)

    workflow.subworkflow(
        func=create_variant_calling_workflow,
        name='create_varcall',
        args=(
            tumour_cells,
            mgd.InputFile('tumour_split.bam', 'region', extensions=['bai']),
            mgd.InputFile('normal_split.bam', 'region', extensions=['bai']),
            mgd.OutputFile(museq_vcf),
            mgd.OutputFile(strelka_snv_vcf),
            mgd.OutputFile(strelka_indel_vcf),
            mgd.OutputFile(snv_h5),
            mgd.OutputFile(meta_yaml),
            config,
            raw_data_dir,
        ),
    )

    return workflow


def variant_calling_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = variant_calling_workflow(args)

    pyp.run(workflow)


def create_variant_calling_workflow(
        tumour_cell_bams,
        tumour_region_bams,
        normal_region_bams,
        museq_vcf,
        strelka_snv_vcf,
        strelka_indel_vcf,
        snv_h5,
        config,
        raw_data_dir,
):
    ctx = {'num_retry': 3,
           'mem_retry_increment': 2,
           'disk_retry_increment': 50,
           'ncpus': 1,
           'docker_image': config['docker']['single_cell_pipeline']}

    basedocker = {'docker_image': config['docker']['single_cell_pipeline']}
    vcftools_docker = {'docker_image': config['docker']['vcftools']}

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.set_filenames('normal_regions.bam', 'region', fnames=normal_region_bams)
    workflow.set_filenames('tumour_cells.bam', 'cell_id', fnames=tumour_cell_bams)
    workflow.set_filenames('tumour_regions.bam', 'region', fnames=tumour_region_bams)

    workflow.transform(
        name="get_regions",
        ctx=dict(mem=config['memory']['low'], **ctx),
        func="single_cell.utils.pysamutils.get_regions_from_reference",
        ret=pypeliner.managed.OutputChunks('region'),
        args=(
            config["ref_genome"],
            config["split_size"],
            config["chromosomes"],
        )
    )

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=list(tumour_cell_bams.keys()),
    )

    workflow.subworkflow(
        name='museq',
        func=mutationseq.create_museq_workflow,
        args=(
            mgd.InputFile('normal_regions.bam', 'region', extensions=['.bai']),
            mgd.InputFile('tumour_regions.bam', 'region', extensions=['.bai']),
            config['ref_genome'],
            mgd.OutputFile(museq_vcf, extensions=['.tbi', '.csi']),
            config,
        ),
    )

    workflow.subworkflow(
        name='strelka',
        func=strelka.create_strelka_workflow,
        args=(
            mgd.InputFile('normal_regions.bam', 'region', extensions=['.bai']),
            mgd.InputFile('tumour_regions.bam', 'region', extensions=['.bai']),
            config['ref_genome'],
            mgd.OutputFile(strelka_indel_vcf, extensions=['.tbi', '.csi']),
            mgd.OutputFile(strelka_snv_vcf, extensions=['.tbi', '.csi']),
            config,
        ),
        kwargs={"chromosomes": config["chromosomes"]}
    )

    workflow.transform(
        name='convert_museq_to_hdf5',
        func="biowrappers.components.io.vcf.tasks.convert_vcf_to_hdf5",
        ctx=dict(mem=2, **ctx),
        args=(
            mgd.InputFile(museq_vcf),
            mgd.TempOutputFile('museq.h5'),
            '/museq/vcf/',
        ),
        kwargs={
            'score_callback': museq_callback,
        }
    )

    workflow.transform(
        name='convert_strelka_to_hdf5',
        func="biowrappers.components.io.vcf.tasks.convert_vcf_to_hdf5",
        ctx=dict(mem=2, **ctx),
        args=(
            mgd.InputFile(strelka_snv_vcf),
            mgd.TempOutputFile('strelka_snv.h5'),
            '/strelka/vcf/',
        ),
        kwargs={
            'score_callback': strelka_snv_callback,
        }
    )

    workflow.transform(
        name='merge_snvs',
        func='biowrappers.components.io.vcf.tasks.merge_vcfs',
        ctx=dict(mem=2, **ctx),
        args=(
            [
                mgd.InputFile(museq_vcf, extensions=['.tbi', '.csi']),
                mgd.InputFile(strelka_snv_vcf, extensions=['.tbi', '.csi']),
            ],
            mgd.TempOutputFile('all.snv.vcf')
        ),
    )

    workflow.transform(
        name='finalise_snvs',
        func="biowrappers.components.io.vcf.tasks.finalise_vcf",
        ctx=dict(mem=2, **ctx),
        args=(
            mgd.TempInputFile('all.snv.vcf'),
            mgd.TempOutputFile('all.snv.vcf.gz', extensions=['.tbi', '.csi'])
        ),
        kwargs={'docker_config': vcftools_docker}
    )

    workflow.subworkflow(
        name='annotate_snvs',
        axes=(),
        ctx=dict(mem=2, **ctx),
        func="biowrappers.pipelines.snv_call_and_annotate.create_annotation_workflow",
        args=(
            config,
            mgd.TempInputFile('all.snv.vcf.gz', extensions=['.tbi', '.csi']),
            mgd.TempOutputFile('snv_annotations.h5'),
            os.path.join(raw_data_dir, 'snv'),
        ),
        kwargs={
            'variant_type': 'snv',
            'docker_config': basedocker,
            'snpeff_docker': vcftools_docker,
            'vcftools_docker': vcftools_docker
        }
    )

    workflow.transform(
        name='build_results_file',
        ctx=dict(mem=config['memory']['high'], **ctx),
        func="biowrappers.components.io.hdf5.tasks.concatenate_tables",
        args=([
                  mgd.TempInputFile('snv_annotations.h5'),
                  mgd.TempInputFile('museq.h5'),
                  mgd.TempInputFile('strelka_snv.h5'),
              ],
              pypeliner.managed.OutputFile(snv_h5),
        ),
        kwargs={
            'drop_duplicates': True,
            'in_memory': False,
        }
    )

    return workflow


def variant_counting_workflow(args):
    config = helpers.load_config(args)

    meta_yaml = os.path.join(args['out_dir'], 'info.yaml')

    bam_files, bai_files = helpers.get_bams(args['input_yaml'])
    vcfs = args['input_vcfs']
    results_file = os.path.join(args['out_dir'], 'results', 'variant_counting', 'counts.h5')

    cellids = helpers.get_samples(args['input_yaml'])

    return create_variant_counting_workflow(vcfs, bam_files, results_file, meta_yaml, config)


def create_variant_counting_workflow(
        vcfs,
        tumour_cell_bams,
        results_h5,
        config,
):
    """ Count variant reads for multiple sets of variants across cells.
    """

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=list(tumour_cell_bams.keys()),
    )

    workflow.transform(
        name='merge_snvs',
        func='biowrappers.components.io.vcf.tasks.merge_vcfs',
        args=(
            [mgd.InputFile(vcf, extensions=['.tbi', '.csi']) for vcf in vcfs],
            mgd.TempOutputFile('all.snv.vcf')
        )
    )

    workflow.transform(
        name='finalise_snvs',
        func="biowrappers.components.io.vcf.tasks.finalise_vcf",
        args=(
            mgd.TempInputFile('all.snv.vcf'),
            mgd.TempOutputFile('all.snv.vcf.gz', extensions=['.tbi'])
        ),
        kwargs={'docker_config': {'docker_image': config['docker']['vcftools']}}
    )

    workflow.subworkflow(
        name='count_alleles',
        func=create_snv_allele_counts_for_vcf_targets_workflow,
        args=(
            mgd.InputFile('tumour_cells.bam', 'cell_id', extensions=['.bai'], fnames=tumour_cell_bams),
            mgd.TempInputFile('all.snv.vcf.gz'),
            mgd.OutputFile(results_h5),
            config['memory'],
        ),
    )

    return workflow


def variant_counting_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = variant_counting_workflow(args)

    pyp.run(workflow)


def variant_calling_multi_sample_workflow(
        config, normal_wgs_bam, tumour_cell_bams, varcall_dir,
        museq_vcf, strelka_snvs, strelka_indels, snv_annotations, snv_counts
):
    keys = [(sample_id, library_id) for (sample_id, library_id, _) in list(tumour_cell_bams.keys())]
    keys = sorted(set(keys))

    museq_vcf = dict([(key, museq_vcf(*key)) for key in keys])
    strelka_snvs = dict([(key, strelka_snvs(*key)) for key in keys])
    strelka_indels = dict([(key, strelka_indels(*key)) for key in keys])
    snv_annotations = dict([(key, snv_annotations(*key)) for key in keys])
    snv_counts = dict([(key, snv_counts(*key)) for key in keys])

    variant_calling_raw_data_template = os.path.join(
        varcall_dir, 'variant_calling_rawdata', '{sample_id}_{library_id}_variant_calling'
    )
    normal_region_bam_template = os.path.join(
        varcall_dir, 'normal_region_bams', 'normal_{region}.bam'
    )
    tumour_region_bam_template = os.path.join(
        varcall_dir, 'tumour_region_bams', '{sample_id}_{library_id}_{region}.bam'
    )

    vcftools_image = {'docker_image': config['variant_calling']['docker']['vcftools']}

    workflow = pypeliner.workflow.Workflow(
        default_ctx={'docker_image': config['multi_sample']['docker']['single_cell_pipeline']}
    )

    workflow.transform(
        name='get_regions',
        ret=mgd.TempOutputObj("get_regions"),
        func=refgenome.get_split_regions,
        args=(
            config["split_bam"]["split_size"],
            config["split_bam"]["ref_genome"]
        )
    )

    workflow.setobj(
        obj=mgd.OutputChunks('region'),
        value=mgd.TempInputObj('get_regions'),
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'library_id', 'cell_id'),
        value=list(tumour_cell_bams.keys()),
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'library_id', 'region'),
        axes=('sample_id', 'library_id',),
        value=mgd.TempInputObj('get_regions'),
    )

    if isinstance(normal_wgs_bam, dict):
        workflow.setobj(
            obj=mgd.OutputChunks('normal_cell_id'),
            value=list(normal_wgs_bam.keys()),
        )
        workflow.set_filenames('normal_cells.bam', 'normal_cell_id', fnames=normal_wgs_bam)
        workflow.subworkflow(
            name="merge_normal_cells",
            func=merge_bams.create_merge_bams_workflow,
            args=(
                mgd.InputFile('normal_cells.bam', 'normal_cell_id', extensions=['.bai']),
                mgd.OutputFile('normal_regions.bam', 'region', axes_origin=[], extensions=['.bai'],
                               template=normal_region_bam_template),
                mgd.TempInputObj('get_regions'),
                config['merge_bams'],
            )
        )
    else:
        workflow.subworkflow(
            name="split_normal",
            func=split_bams.create_split_workflow,
            args=(
                mgd.InputFile(normal_wgs_bam, extensions=['.bai']),
                mgd.OutputFile('normal_regions.bam', 'region', extensions=['.bai'], axes_origin=[],
                               template=normal_region_bam_template),
                pypeliner.managed.TempInputObj('region'),
                config['split_bam'],
            ),
            kwargs={"by_reads": False}
        )

    workflow.subworkflow(
        name="split_merge_tumour",
        axes=('sample_id', 'library_id',),
        func=merge_bams.create_merge_bams_workflow,
        args=(
            mgd.InputFile('tumour_all_cells.bam', 'sample_id', 'library_id', 'cell_id', fnames=tumour_cell_bams,
                          extensions=['.bai'], axes_origin=[]),
            mgd.OutputFile('tumour_regions.bam', 'sample_id', 'library_id', 'region', axes_origin=[],
                           extensions=['.bai'], template=tumour_region_bam_template),
            mgd.TempInputObj('get_regions'),
            config['merge_bams'],
        )
    )

    workflow.subworkflow(
        name='variant_calling',
        func=create_variant_calling_workflow,
        axes=('sample_id', 'library_id',),
        args=(
            mgd.InputFile('tumour_all_cells.bam', 'sample_id', 'library_id', 'cell_id', extensions=['.bai'],
                          fnames=tumour_cell_bams),
            mgd.InputFile('tumour_regions.bam', 'sample_id', 'library_id', 'region', extensions=['.bai'],
                          template=tumour_region_bam_template),
            mgd.InputFile('normal_regions.bam', 'region', extensions=['.bai'], template=normal_region_bam_template),
            mgd.OutputFile('museq.vcf', 'sample_id', 'library_id', extensions=['.tbi', '.csi'], fnames=museq_vcf),
            mgd.OutputFile('strelka_snv.vcf', 'sample_id', 'library_id', extensions=['.tbi', '.csi'],
                           fnames=strelka_snvs),
            mgd.OutputFile('strelka_indel.vcf', 'sample_id', 'library_id', extensions=['.tbi', '.csi'],
                           fnames=strelka_indels),
            mgd.OutputFile('snv_annotations.h5', 'sample_id', 'library_id', fnames=snv_annotations),
            config['variant_calling'],
            mgd.Template(variant_calling_raw_data_template, 'sample_id', 'library_id'),
        ),
    )

    workflow.transform(
        name='merge_museq_snvs',
        func='biowrappers.components.io.vcf.tasks.concatenate_vcf',
        args=(
            mgd.InputFile('museq.vcf', 'sample_id', 'library_id', axes_origin=[], extensions=['.tbi', '.csi'],
                          fnames=museq_vcf),
            mgd.TempOutputFile('museq.vcf.gz', extensions=['.tbi', '.csi']),
        ),
        kwargs={
            'allow_overlap': True,
            'docker_config': vcftools_image,
        },
    )

    workflow.transform(
        name='merge_strelka_snvs',
        func='biowrappers.components.io.vcf.tasks.concatenate_vcf',
        args=(
            mgd.InputFile('strelka_snv.vcf', 'sample_id', 'library_id',
                          axes_origin=[], extensions=['.tbi', '.csi'], fnames=strelka_snvs),
            mgd.TempOutputFile('strelka_snv.vcf.gz', extensions=['.tbi', '.csi']),
        ),
        kwargs={
            'allow_overlap': True,
            'docker_config': vcftools_image,
        },
    )

    workflow.subworkflow(
        name='variant_counting',
        func=create_variant_counting_workflow,
        axes=('sample_id', 'library_id',),
        args=(
            [
                mgd.TempInputFile('museq.vcf.gz', extensions=['.tbi', '.csi']),
                mgd.TempInputFile('strelka_snv.vcf.gz', extensions=['.tbi', '.csi']),
            ],
            mgd.InputFile('tumour_all_cells.bam', 'sample_id', 'library_id', 'cell_id', extensions=['.bai'],
                          fnames=tumour_cell_bams),
            mgd.OutputFile('snv_counts.h5', 'sample_id', 'library_id', fnames=snv_counts),
            config['variant_calling'],
        ),
    )

    return workflow
