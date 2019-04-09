'''
Created on Feb 22, 2018

@author: dgrewal
'''

import os
import pypeliner
import pypeliner.managed as mgd
from workflows import mutationseq
from workflows import strelka
import single_cell
from single_cell.utils import helpers

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

    ctx={'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1,
         'mem': config["memory"]['low'], 'docker_image': baseimage}
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    if isinstance(normal_bams, dict) and isinstance(tumour_bams, dict):
        assert list(normal_bams.keys()) == list(tumour_bams.keys()), 'keys for tumour and normal bams should be the same'
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


def create_variant_calling_workflow(
        tumour_cell_bams,
        tumour_region_bams,
        normal_region_bams,
        museq_vcf,
        strelka_snv_vcf,
        strelka_indel_vcf,
        snv_h5,
        meta_yaml,
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

    workflow.subworkflow(
        name='count_alleles',
        func=create_snv_allele_counts_for_vcf_targets_workflow,
        ctx=ctx,
        args=(
            mgd.InputFile('tumour_cells.bam', 'cell_id', extensions=['.bai'], fnames=tumour_cell_bams),
            mgd.TempInputFile('all.snv.vcf.gz'),
            mgd.TempOutputFile('snv_counts.h5'),
            config['memory'],
        ),
    )

    workflow.transform(
        name='build_results_file',
        ctx=dict(mem=config['memory']['high'], **ctx),
        func="biowrappers.components.io.hdf5.tasks.concatenate_tables",
        args=([
                  mgd.TempInputFile('snv_counts.h5'),
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

    normals = {k: helpers.format_file_yaml(v) for k,v in normal_region_bams.iteritems()}
    tumours = {k: helpers.format_file_yaml(v) for k,v in tumour_region_bams.iteritems()}
    cells = {k: helpers.format_file_yaml(v) for k,v in tumour_cell_bams.iteritems()}
    inputs = {'normal': normals, 'tumour': tumours, 'cells':cells}

    metadata = {
        'variant_calling': {
            'name': 'variant_calling',
            'version': single_cell.__version__,
            'containers': config['docker'],
            'output_datasets': None,
            'input_datasets': inputs,
            'results': {'variant_calling_data': helpers.format_file_yaml(snv_h5)}
        }
    }

    # TODO: will download results unnecessarily on cloud
    workflow.transform(
        name='generate_meta_yaml',
        ctx=dict(mem=config['memory']['med']),
        func="single_cell.utils.helpers.write_to_yaml",
        args=(
            mgd.OutputFile(meta_yaml),
            metadata
        )
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
        meta_yaml,
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

    # TODO: will download results unnecessarily on cloud
    workflow.transform(
        name='generate_meta_yaml',
        func="single_cell.utils.helpers.write_to_yaml",
        args=(
            mgd.OutputFile(meta_yaml),
            {
                'name': 'variant_counting',
                'version': single_cell.__version__,
                'output_datasets': None,
                'input_datasets': {
                    'tumour_cell_bam': tumour_cell_bams,
                },
                'results': {
                    'snv_counts': mgd.InputFile(results_h5),
                },
            },
        )
    )

    return workflow
