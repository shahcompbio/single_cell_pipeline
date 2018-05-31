'''
Created on Feb 22, 2018

@author: dgrewal
'''

import os
import pypeliner
import pypeliner.managed as mgd
from workflows import germline
from workflows import split_bams
from single_cell.utils import helpers
from single_cell.workflows.germline import tasks
import variant_calling


def germline_calling_workflow(workflow, args):

    config = helpers.load_config(args)

    bam_files, bai_files  = helpers.get_bams(args['input_yaml'])

    sampleids = helpers.get_samples(args['input_yaml'])

    normal_bam_template = args["input_template"]
    normal_bai_template = args["input_template"] + ".bai"

    if "{reads}" in normal_bam_template:
        raise ValueError("input template for germline calling only support region based splits")


    varcalls_dir = os.path.join(
        args['out_dir'], 'results', 'germline_calling')

    samtools_germline_vcf = os.path.join(varcalls_dir, 'raw', 'samtools_germline.vcf.gz')
    snpeff_vcf_filename = os.path.join(varcalls_dir, 'snpeff.vcf')
    normal_genotype_filename = os.path.join(varcalls_dir, 'raw', 'normal_genotype.h5')
    mappability_filename = os.path.join(varcalls_dir, 'raw', 'mappability.h5')
    counts_template = os.path.join(varcalls_dir, 'counts', 'raw', 'counts.h5')
    germline_h5_filename = os.path.join(varcalls_dir, 'germline.h5')

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=bam_files.keys(),
    )
 
    workflow.transform(
        name="get_regions",
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 },
        func=helpers.get_regions_from_reference,
        ret=pypeliner.managed.OutputChunks('region'),
        args=(
            config["ref_genome"],
            config["split_size"],
            config["chromosomes"],
        )
    )
 
    workflow.subworkflow(
        name='samtools_germline',
        func=germline.create_samtools_germline_workflow,
        args=(
            mgd.InputFile("normal.split.bam", "region", template=normal_bam_template),
            mgd.InputFile("normal.split.bam.bai", "region", template=normal_bai_template),
            config['ref_genome'],
            mgd.OutputFile(samtools_germline_vcf, extensions=['.tbi']),
            config,
        ),
    )

    workflow.subworkflow(
        name='annotate_mappability',
        func="biowrappers.components.variant_calling.mappability.create_vcf_mappability_annotation_workflow",
        args=(
            config['databases']['mappability']['local_path'],
            mgd.InputFile(samtools_germline_vcf, extensions=['.tbi']),
            mgd.OutputFile(mappability_filename),
        ),
    )

    workflow.transform(
        name='annotate_genotype',
        func=tasks.annotate_normal_genotype,
        args=(
            mgd.InputFile(samtools_germline_vcf, extensions=['.tbi']),
            mgd.OutputFile(normal_genotype_filename),
            config["chromosomes"],
        ),
    )

    workflow.subworkflow(
        name='snpeff',
        func=biowrappers.components.variant_calling.snpeff.create_snpeff_annotation_workflow,
        args=(
            config['databases']['snpeff']['db'],
            mgd.InputFile(samtools_germline_vcf, extensions=['.tbi']),
            mgd.OutputFile(snpeff_vcf_filename),
        ),
        kwargs={
            'hdf5_output': False,
        }
    )

    workflow.subworkflow(
        name='read_counts',
        func=variant_calling.create_snv_allele_counts_for_vcf_targets_workflow,
        args=(
            config,
            mgd.InputFile('tumour.bam', 'cell_id', fnames=bam_files),
            mgd.InputFile('tumour.bam.bai', 'cell_id', fnames=bai_files),
            mgd.InputFile(samtools_germline_vcf, extensions=['.tbi']),
            mgd.OutputFile(counts_template),
        ),
        kwargs={
            'table_name': '/germline_allele_counts',
        },
    )

    workflow.transform(
        name='build_results_file',
        func=biowrappers.components.io.hdf5.tasks.concatenate_tables,
        args=([
                mgd.InputFile(counts_template),
                mgd.InputFile(mappability_filename),
                mgd.InputFile(normal_genotype_filename),
            ],
            pypeliner.managed.OutputFile(germline_h5_filename),
        ),
        kwargs={
            'drop_duplicates' : True,
        }
    )

    return workflow

