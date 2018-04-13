'''
Created on Feb 22, 2018

@author: dgrewal
'''

import os
import pypeliner
import pypeliner.managed as mgd
import biowrappers.components.io.vcf.tasks
import biowrappers.components.io.hdf5.tasks
import biowrappers.pipelines.snv_call_and_annotate
from workflows import snv_postprocessing
from workflows import mutationseq 
from workflows import split_bams
from workflows import strelka
from single_cell.utils import helpers
from single_cell.utils import vcfutils


default_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']


def create_snv_allele_counts_for_vcf_targets_workflow(
        config,
        bam_files,
        bai_files,
        vcf_file,
        out_file,
        chromosomes=default_chromosomes,
        count_duplicates=False,
        min_bqual=0,
        min_mqual=0,
        split_size=int(1e7),
        table_name='snv_allele_counts',
        vcf_to_bam_chrom_map=None):

    workflow = pypeliner.workflow.Workflow(default_ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 })

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=bam_files.keys(),
    )

    workflow.transform(
        name='get_snv_allele_counts_for_vcf_targets',
        axes=('cell_id',),
        func=biowrappers.components.variant_calling.snv_allele_counts.tasks.get_snv_allele_counts_for_vcf_targets,
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
        }
    )

    workflow.transform(
        name='merge_snv_allele_counts',
        func=biowrappers.components.io.hdf5.tasks.concatenate_tables,
        args=(
            mgd.TempInputFile('counts.h5', 'cell_id'),
            mgd.OutputFile(out_file),
        ),
    )

    return workflow


def variant_calling_workflow(workflow, args):

    args = helpers.generate_configs_in_temp(args)

    config = helpers.load_config(args)

    bam_files, bai_files  = helpers.get_bams(args['input_yaml'])

    cellids = helpers.get_samples(args['input_yaml'])

    varcalls_dir = os.path.join(args['out_dir'], 'results',
                                'variant_calling')

    museq_vcf = os.path.join(varcalls_dir, 'museq_snv.vcf')
    strelka_snv_vcf = os.path.join(varcalls_dir, 'strelka_snv.vcf')
    strelka_indel_vcf = os.path.join(varcalls_dir, 'strelka_indel.vcf')
    snv_h5_filename = os.path.join(varcalls_dir, 'snv_annotations.h5')

    wgs_bam_template = args["tumour_template"]
    wgs_bai_template = args["tumour_template"] + ".bai"

    normal_bam_template = args["normal_template"]
    normal_bai_template = args["normal_template"] + ".bai"

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cellids,
    )

    workflow.transform(
        name="get_regions",
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 },
        func=helpers.get_bam_regions,
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
            mgd.InputFile("normal.split.bam", "region", template=normal_bam_template, axes_origin=[]),
            mgd.InputFile("normal.split.bam.bai", "region", template=normal_bai_template, axes_origin=[]),
            mgd.InputFile("merged_bam", "region", template=wgs_bam_template, axes_origin=[]),
            mgd.InputFile("merged_bam", "region", template=wgs_bai_template, axes_origin=[]),
            config['ref_genome'],
            mgd.OutputFile(museq_vcf),
            config,
        ),
    )

    workflow.subworkflow(
        name='strelka',
        func=strelka.create_strelka_workflow,
        args=(
            mgd.InputFile("normal.split.bam", "region", template=normal_bam_template, axes_origin=[]),
            mgd.InputFile("normal.split.bam.bai", "region", template=normal_bai_template, axes_origin=[]),
            mgd.InputFile("merged_bam", "region", template=wgs_bam_template, axes_origin=[]),
            mgd.InputFile("merged_bam", "region", template=wgs_bai_template, axes_origin=[]),
            config['ref_genome'],
            mgd.OutputFile(strelka_indel_vcf),
            mgd.OutputFile(strelka_snv_vcf),
            config,
        ),
    )

    workflow.transform(
        name='merge_snvs',
        func=vcfutils.merge_vcfs,
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
        func=biowrappers.components.io.vcf.tasks.finalise_vcf,
        args=(
            mgd.TempInputFile('all.snv.vcf'),
            mgd.TempOutputFile('all.snv.vcf.gz')
        )
    )

    workflow.subworkflow(
        name='annotate_snvs',
        axes=(),
        func=biowrappers.pipelines.snv_call_and_annotate.create_annotation_workflow,
        args=(
            config,
            mgd.TempInputFile('all.snv.vcf.gz'),
            mgd.TempOutputFile('snv_annotations.h5'),
            os.path.join(varcalls_dir, 'snv'),
        ),
        kwargs={
            'variant_type': 'snv'
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
        )
    )

    workflow.transform(
        name='build_results_file',
        func=biowrappers.components.io.hdf5.tasks.concatenate_tables,
        args=([
                mgd.TempInputFile('snv_counts.h5'),
                mgd.TempInputFile('snv_annotations.h5'),
            ],
            pypeliner.managed.OutputFile(snv_h5_filename),
        ),
        kwargs={
            'drop_duplicates' : True,
        }
    )

    return workflow

