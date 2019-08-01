import pypeliner

default_chromosomes = [str(a) for a in range(1, 23)] + ['X', 'Y']


def create_vardict_paired_sample_workflow(
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        out_file,
        config,
        chromosomes=default_chromosomes,
        java=False,
        min_allele_frequency=0.01,
        remove_duplicate_reads=False,
        sample_names=None):
    workflow = pypeliner.workflow.Workflow(
        default_ctx={
            'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2, 'ncpus': 1, 'disk_retry_increment': 50
        })

    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('region'),
        value=list(normal_bam_file.keys()),
    )

    workflow.transform(
        name='run_vardict',
        axes=('region',),
        ctx={'mem': 12},
        func="biowrappers.components.variant_calling.vardict.tasks.run_paired_sample_vardict",
        args=(
            pypeliner.managed.InputFile('normal.bam', 'region', fnames=normal_bam_file, extensions=['.bai']),
            pypeliner.managed.InputFile('tumour.bam', 'region', fnames=tumour_bam_file, extensions=['.bai']),
            ref_genome_fasta_file,
            pypeliner.managed.InputInstance('region'),
            pypeliner.managed.TempOutputFile('result.vcf', 'region'),
        ),
        kwargs={
            'java': java,
            'min_allele_frequency': min_allele_frequency,
            'remove_duplicate_reads': remove_duplicate_reads,
            'sample_names': sample_names,
        },
    )

    workflow.transform(
        name='compress_tmp',
        axes=('region',),
        func="biowrappers.components.io.vcf.tasks.compress_vcf",
        args=(
            pypeliner.managed.TempInputFile('result.vcf', 'region'),
            pypeliner.managed.TempOutputFile('result.vcf.gz', 'region'),
        ),
    )

    workflow.transform(
        name='concatenate_vcf',
        func="biowrappers.components.io.vcf.tasks.concatenate_vcf",
        args=(
            pypeliner.managed.TempInputFile('result.vcf.gz', 'region'),
            pypeliner.managed.OutputFile(out_file),
        ),
    )

    return workflow
