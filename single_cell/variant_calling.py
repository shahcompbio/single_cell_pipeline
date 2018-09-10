'''
Created on Feb 22, 2018

@author: dgrewal
'''

import os
import pypeliner
import pypeliner.managed as mgd
from workflows import mutationseq 
from workflows import strelka
from single_cell.utils import helpers


default_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']


def create_snv_allele_counts_for_vcf_targets_workflow(
        config,
        bam_files,
        bai_files,
        vcf_file,
        out_file,
        docker_config=None,
        chromosomes=default_chromosomes,
        count_duplicates=False,
        min_bqual=0,
        min_mqual=0,
        split_size=int(1e7),
        table_name='snv_allele_counts',
        vcf_to_bam_chrom_map=None):

    ctx = {'mem': 2, 'num_retry': 3,
           'mem_retry_increment': 2,
           'pool_id': config['pools']['standard'], 'ncpus': 1}
    if docker_config:
        ctx.update(docker_config)

    workflow = pypeliner.workflow.Workflow(default_ctx=ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=bam_files.keys(),
    )

    workflow.transform(
        name='get_snv_allele_counts_for_vcf_targets',
        axes=('cell_id',),
        func="biowrappers.components.variant_calling.snv_allele_counts.tasks.get_snv_allele_counts_for_vcf_targets",
        args=(
            mgd.InputFile('tumour.bam', 'cell_id', fnames=bam_files),
            mgd.InputFile('tumour.bam.bai', 'cell_id', fnames=bai_files),
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
        ctx={'mem': config["memory"]['high'], 'pool_id': config['pools']['highmem'], 'ncpus':1},
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


def variant_calling_workflow(workflow, args):

    config = helpers.load_config(args)

    ctx = {'num_retry': 3,
           'mem_retry_increment': 2,
           'ncpus': 1}
    docker_ctx = helpers.build_docker_args(config['docker'], 'single_cell_pipeline')
    ctx.update(docker_ctx)

    bam_files, bai_files = helpers.get_bams(args['input_yaml'])

    cellids = helpers.get_samples(args['input_yaml'])

    varcalls_dir = os.path.join(args['out_dir'], 'results',
                                'variant_calling')

    museq_vcf = os.path.join(varcalls_dir, 'museq_snv.vcf.gz')
    strelka_snv_vcf = os.path.join(varcalls_dir, 'strelka_snv.vcf.gz')
    strelka_indel_vcf = os.path.join(varcalls_dir, 'strelka_indel.vcf.gz')
    snv_h5_filename = os.path.join(varcalls_dir, 'snv_annotations.h5')

    wgs_bam_template = args["tumour_template"]
    wgs_bai_template = args["tumour_template"] + ".bai"

    normal_bam_template = args["normal_template"]
    normal_bai_template = args["normal_template"] + ".bai"

    singlecellimage = config['docker']['images']['single_cell_pipeline']

    if "{reads}" in normal_bam_template or "{reads}" in wgs_bam_template:
        raise ValueError("input template for variant calling only supports region based splits")


    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cellids,
    )

    workflow.transform(
        name="get_regions",
        ctx=dict(mem=2, pool_id=config['pools']['standard'], **ctx),
        func="single_cell.utils.pysamutils.get_regions_from_reference",
        ret=mgd.OutputChunks('region'),
        args=(
              config["ref_genome"],
              config["split_size"],
              config["chromosomes"],
        )
    )

    workflow.subworkflow(
        name='museq',
        func=mutationseq.create_museq_workflow,
        args=(
            mgd.InputFile("normal.split.bam", "region", template=normal_bam_template),
            mgd.InputFile("normal.split.bam.bai", "region", template=normal_bai_template),
            mgd.InputFile("merged_bam", "region", template=wgs_bam_template),
            mgd.InputFile("merged_bam", "region", template=wgs_bai_template),
            config['ref_genome'],
            mgd.OutputFile(museq_vcf),
            config,
        ),
    )

    workflow.subworkflow(
        name='strelka',
        func=strelka.create_strelka_workflow,
        args=(
            mgd.InputFile("normal.split.bam", "region", template=normal_bam_template),
            mgd.InputFile("normal.split.bam.bai", "region", template=normal_bai_template),
            mgd.InputFile("merged_bam", "region", template=wgs_bam_template),
            mgd.InputFile("merged_bam", "region", template=wgs_bai_template),
            config['ref_genome'],
            mgd.OutputFile(strelka_indel_vcf),
            mgd.OutputFile(strelka_snv_vcf),
            config,
        ),
        kwargs={"chromosomes":config["chromosomes"]}
    )

    workflow.transform(
        name='convert_museq_to_hdf5',
        func="biowrappers.components.io.vcf.tasks.convert_vcf_to_hdf5",
        ctx=dict(mem=2, pool_id=config['pools']['standard'], **ctx),
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
        ctx=dict(mem=2, pool_id=config['pools']['standard'], **ctx),
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
        func="single_cell.utils.vcfutils.merge_vcfs",
        ctx=dict(mem=2, pool_id=config['pools']['standard'], **ctx),
        args=(
            [
                mgd.InputFile(museq_vcf),
                mgd.InputFile(strelka_snv_vcf),
            ],
            mgd.TempOutputFile('all.snv.vcf')
        )
    )

    workflow.transform(
        name='finalise_snvs',
        func="biowrappers.components.io.vcf.tasks.finalise_vcf",
        ctx=dict(mem=2, pool_id=config['pools']['standard'], **ctx),
        args=(
            mgd.TempInputFile('all.snv.vcf'),
            mgd.TempOutputFile('all.snv.vcf.gz', extensions=['.tbi'])
        )
    )

    workflow.subworkflow(
        name='annotate_snvs',
        axes=(),
        func="biowrappers.pipelines.snv_call_and_annotate.create_annotation_workflow",
        args=(
            config,
            mgd.TempInputFile('all.snv.vcf.gz'),
            mgd.TempOutputFile('snv_annotations.h5'),
            os.path.join(varcalls_dir, 'snv'),
        ),
        kwargs={
            'variant_type': 'snv',
            'docker_config': helpers.build_docker_args(config['docker'], 'single_cell_pipeline')
        }
    )

    workflow.subworkflow(
        name='count_alleles',
        func=create_snv_allele_counts_for_vcf_targets_workflow,
        args=(
            config,
            mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_files),
            mgd.InputFile('bam_markdups_index', 'cell_id', fnames=bai_files),
            mgd.TempInputFile('all.snv.vcf.gz'),
            mgd.TempOutputFile('snv_counts.h5'),
        ),
        kwargs={'chromosomes': config['chromosomes'],
                'docker_config': config['docker']}
    )

    workflow.transform(
        name='build_results_file',
        ctx=dict(mem=config['memory']['high'], pool_id=config['pools']['highmem'], **ctx),
        func="biowrappers.components.io.hdf5.tasks.concatenate_tables",
        args=([
                mgd.TempInputFile('snv_counts.h5'),
                mgd.TempInputFile('snv_annotations.h5'),
                mgd.TempInputFile('museq.h5'),
                mgd.TempInputFile('strelka_snv.h5'),
            ],
            pypeliner.managed.OutputFile(snv_h5_filename),
        ),
        kwargs={
            'drop_duplicates' : True,
            'in_memory' : False,
        }
    )

    return workflow


def variant_counting_workflow(workflow, args):

    config = helpers.load_config(args)

    bam_files, bai_files  = helpers.get_bams(args['input_yaml'])
    vcfs = args['input_vcfs']
    results_file = os.path.join(args['out_dir'], 'results', 'variant_counting', 'counts.h5')

    cellids = helpers.get_samples(args['input_yaml'])

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cellids,
    )

    workflow.transform(
        name='merge_snvs',
        func="single_cell.utils.vcfutils.merge_vcfs",
        args=(
            [mgd.InputFile(vcf) for vcf in vcfs],
            mgd.TempOutputFile('all.snv.vcf')
        )
    )

    workflow.transform(
        name='finalise_snvs',
        func="biowrappers.components.io.vcf.tasks.finalise_vcf",
        args=(
            mgd.TempInputFile('all.snv.vcf'),
            mgd.TempOutputFile('all.snv.vcf.gz', extensions=['.tbi'])
        )
    )

    workflow.subworkflow(
        name='count_alleles',
        func=create_snv_allele_counts_for_vcf_targets_workflow,
        args=(
            config,
            mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_files),
            mgd.InputFile('bam_markdups_index', 'cell_id', fnames=bai_files),
            mgd.TempInputFile('all.snv.vcf.gz'),
            mgd.OutputFile(results_file),
        ),
        kwargs={
            'docker_config': helpers.build_docker_args(config['docker'], 'single_cell_pipeline')
        },
    )

    return workflow

