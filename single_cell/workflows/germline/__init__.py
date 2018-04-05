from pypeliner.workflow import Workflow

import csv
import pypeliner
import pysam
import single_cell.workflows.strelka.vcf_tasks as vcf_tasks
import tasks
import os


default_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']


def create_samtools_germline_workflow(
        normal_bam_files,
        normal_bai_files,
        ref_genome_fasta_file,
        vcf_file,
        parsed_csv,
        config,
        regions,
        chromosomes=default_chromosomes,
        split_size=int(1e7)):

    workflow = Workflow()

    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('regions'),
        value=regions,
    )

    workflow.transform(
        name='run_samtools_variant_calling',
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2, 'ncpus':1, 'pool_id': config['pools']['standard']},
        axes=('regions',),
        func=tasks.run_samtools_variant_calling,
        args=(
            pypeliner.managed.InputFile('normal.bam', 'regions', fnames=normal_bam_files),
            pypeliner.managed.InputFile('normal.bam.bai', 'regions', fnames=normal_bai_files),
            pypeliner.managed.InputFile(ref_genome_fasta_file),
            pypeliner.managed.TempOutputFile('variants.vcf.gz', 'regions'),
        ),
        kwargs={
            'region': pypeliner.managed.InputInstance('regions'),
        },
    )

    workflow.transform(
        name='concatenate_variants',
        ctx={'mem': 2},
        func=vcf_tasks.concatenate_vcf,
        args=(
            pypeliner.managed.TempInputFile('variants.vcf.gz', 'regions'),
            pypeliner.managed.OutputFile(vcf_file),
        ),
    )

    workflow.transform(
        name='parse_samtools_vcf',
        ctx={'mem': config["memory"]['low'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.parse_samtools_vcf,
        args=(
            mgd.InputFile(vcf_file),
            mgd.OutputFile(parsed_csv),
        ),
    )

