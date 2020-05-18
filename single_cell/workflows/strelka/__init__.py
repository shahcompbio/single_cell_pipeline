from pypeliner.workflow import Workflow

import pypeliner
from single_cell.workflows.strelka.dtypes import dtypes


default_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']


def strelka_snv_callback(record):
    return record.INFO['QSS']


def create_strelka_workflow(
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        indel_vcf_file,
        snv_vcf_file,
        snv_csv_file,
        config,
        chromosomes=default_chromosomes,
        split_size=int(1e7),
        use_depth_thresholds=True):
    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1,
           'num_retry': 3, 'docker_image': config['docker']['single_cell_pipeline']}

    strelka_docker = {'docker_image': config['docker']['strelka']}
    vcftools_docker = {'docker_image': config['docker']['vcftools']}

    regions = list(normal_bam_file.keys())
    assert set(tumour_bam_file.keys()) == set(regions)

    workflow = Workflow(ctx=ctx)

    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('chrom'),
        value=chromosomes,
    )

    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('region'),
        value=regions,
    )

    workflow.transform(
        name='count_fasta_bases',
        ctx=dict(mem=2),
        func="single_cell.workflows.strelka.tasks.count_fasta_bases",
        args=(
            ref_genome_fasta_file,
            pypeliner.managed.TempOutputFile('ref_base_counts.tsv'),
            strelka_docker
        )
    )

    workflow.transform(
        name="get_chrom_sizes",
        ctx=dict(mem=2),
        func="single_cell.workflows.strelka.tasks.get_known_chromosome_sizes",
        ret=pypeliner.managed.TempOutputObj('known_sizes'),
        args=(
            pypeliner.managed.TempInputFile('ref_base_counts.tsv'),
            chromosomes
        )
    )

    workflow.transform(
        name='call_somatic_variants',
        ctx=dict(mem=4, disk=40),
        func="single_cell.workflows.strelka.tasks.call_somatic_variants",
        axes=('region',),
        args=(
            pypeliner.managed.InputFile("normal.split.bam", "region", fnames=normal_bam_file, extensions=['.bai']),
            pypeliner.managed.InputFile("merged_bam", "region", fnames=tumour_bam_file, extensions=['.bai']),
            pypeliner.managed.TempInputObj('known_sizes'),
            ref_genome_fasta_file,
            pypeliner.managed.TempOutputFile('somatic.indels.unfiltered.vcf', 'region'),
            pypeliner.managed.TempOutputFile('somatic.indels.unfiltered.vcf.window', 'region'),
            pypeliner.managed.TempOutputFile('somatic.snvs.unfiltered.vcf', 'region'),
            pypeliner.managed.TempOutputFile('strelka.stats', 'region'),
            pypeliner.managed.InputInstance("region"),
            strelka_docker
        ),
    )

    workflow.transform(
        name='add_indel_filters',
        axes=('chrom',),
        ctx=dict(mem=4),
        func="single_cell.workflows.strelka.tasks.filter_indel_file_list",
        args=(
            pypeliner.managed.TempInputFile('somatic.indels.unfiltered.vcf', 'region'),
            pypeliner.managed.TempInputFile('strelka.stats', 'region'),
            pypeliner.managed.TempInputFile('somatic.indels.unfiltered.vcf.window', 'region'),
            pypeliner.managed.TempOutputFile('somatic.indels.filtered.vcf', 'chrom'),
            pypeliner.managed.InputInstance("chrom"),
            pypeliner.managed.TempInputObj('known_sizes'),
            regions
        ),
        kwargs={'use_depth_filter': use_depth_thresholds}
    )

    workflow.transform(
        name='add_snv_filters',
        axes=('chrom',),
        ctx=dict(mem=4),
        func="single_cell.workflows.strelka.tasks.filter_snv_file_list",
        args=(
            pypeliner.managed.TempInputFile('somatic.snvs.unfiltered.vcf', 'region'),
            pypeliner.managed.TempInputFile('strelka.stats', 'region'),
            pypeliner.managed.TempOutputFile('somatic.snvs.filtered.vcf', 'chrom'),
            pypeliner.managed.InputInstance("chrom"),
            pypeliner.managed.TempInputObj('known_sizes'),
            regions,
        ),
        kwargs={'use_depth_filter': use_depth_thresholds}
    )

    workflow.transform(
        name='merge_indels',
        ctx=dict(mem=4),
        func="single_cell.workflows.strelka.vcf_tasks.concatenate_vcf",
        args=(
            pypeliner.managed.TempInputFile('somatic.indels.filtered.vcf', 'chrom'),
            pypeliner.managed.TempOutputFile('somatic.indels.filtered.vcf.gz'),
            pypeliner.managed.TempSpace("merge_indels_temp"),
            vcftools_docker
        )
    )

    workflow.transform(
        name='merge_snvs',
        ctx=dict(mem=4),
        func="single_cell.workflows.strelka.vcf_tasks.concatenate_vcf",
        args=(
            pypeliner.managed.TempInputFile('somatic.snvs.filtered.vcf', 'chrom'),
            pypeliner.managed.TempOutputFile('somatic.snvs.filtered.vcf.gz'),
            pypeliner.managed.TempSpace("merge_snvs_temp"),
            vcftools_docker
        )
    )

    workflow.transform(
        name='filter_indels',
        ctx=dict(mem=4),
        func="single_cell.workflows.strelka.vcf_tasks.filter_vcf",
        args=(
            pypeliner.managed.TempInputFile('somatic.indels.filtered.vcf.gz'),
            pypeliner.managed.TempOutputFile('somatic.indels.passed.vcf')
        )
    )

    workflow.transform(
        name='filter_snvs',
        ctx=dict(mem=4),
        func="single_cell.workflows.strelka.vcf_tasks.filter_vcf",
        args=(
            pypeliner.managed.TempInputFile('somatic.snvs.filtered.vcf.gz'),
            pypeliner.managed.TempOutputFile('somatic.snvs.passed.vcf')
        )
    )

    workflow.transform(
        name='finalise_indels',
        ctx=dict(mem=4),
        func="single_cell.workflows.strelka.vcf_tasks.finalise_vcf",
        args=(
            pypeliner.managed.TempInputFile('somatic.indels.passed.vcf'),
            pypeliner.managed.OutputFile(indel_vcf_file, extensions=['.tbi', '.csi']),
            vcftools_docker
        )
    )

    workflow.transform(
        name='finalise_snvs',
        ctx=dict(mem=2),
        func="single_cell.workflows.strelka.vcf_tasks.finalise_vcf",
        args=(
            pypeliner.managed.TempInputFile('somatic.snvs.passed.vcf'),
            pypeliner.managed.OutputFile(snv_vcf_file, extensions=['.tbi', '.csi']),
            vcftools_docker
        )
    )

    workflow.transform(
        name='convert_strelka_to_csv',
        func="biowrappers.components.io.vcf.tasks.convert_vcf_to_csv",
        ctx=ctx,
        args=(
            pypeliner.managed.InputFile(snv_vcf_file),
            pypeliner.managed.TempOutputFile('strelka_snv.csv'),
        ),
        kwargs={
            'score_callback': strelka_snv_callback,
        }
    )

    workflow.transform(
        name='prep_strelka_csv',
        func='single_cell.utils.csvutils.rewrite_csv_file',
        args=(
            pypeliner.managed.TempInputFile('strelka_snv.csv'),
            pypeliner.managed.OutputFile(snv_csv_file, extensions=['.yaml'])
        ),
        kwargs={'dtypes': dtypes()['snv_strelka']}
    )


    return workflow
