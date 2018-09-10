from pypeliner.workflow import Workflow
import pypeliner


default_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']


def create_samtools_germline_workflow(
        normal_bam_files,
        normal_bai_files,
        ref_genome_fasta_file,
        vcf_file,
        config,
        chromosomes=default_chromosomes,
        base_docker=None,
        samtools_docker=None,
        vcftools_docker=None
):

    ctx = {'mem': config["memory"]['low'],
           'pool_id': config['pools']['standard'],
           'mem_retry_increment': 2,
           'ncpus': 1}
    if base_docker:
        ctx.update(base_docker)

    regions = normal_bam_files.keys()

    workflow = Workflow()

    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('regions'),
        value=regions,
    )

    workflow.transform(
        name='run_samtools_variant_calling',
        ctx=ctx,
        axes=('regions',),
        func="single_cell.workflows.germline.tasks.run_samtools_variant_calling",
        args=(
            pypeliner.managed.InputFile('normal.split.bam', 'regions', fnames=normal_bam_files),
            pypeliner.managed.InputFile('normal.split.bam.bai', 'regions', fnames=normal_bai_files),
            ref_genome_fasta_file,
            pypeliner.managed.TempOutputFile('variants.vcf.gz', 'regions'),
        ),
        kwargs={
            'region': pypeliner.managed.InputInstance('regions'),
            'samtools_docker': samtools_docker,
            'vcftools_docker': samtools_docker
        },
    )
  
    workflow.transform(
        name='concatenate_variants',
        ctx=ctx,
        func="single_cell.workflows.strelka.vcf_tasks.concatenate_vcf",
        args=(
            pypeliner.managed.TempInputFile('variants.vcf.gz', 'regions'),
            pypeliner.managed.OutputFile(vcf_file, extensions=['.tbi']),
            pypeliner.managed.TempSpace("merge_variants_germline"),
        ),
        kwargs={'docker_config': vcftools_docker}
    )

    return workflow
