'''
Created on Jul 24, 2017

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers
def create_museq_workflow(
        normal_bam, normal_bai, tumour_bam, tumour_bai, ref_genome, snv_vcf,
        config):

    ctx = {'mem_retry_increment': 2, 'ncpus': 1}
    docker_ctx = helpers.get_container_ctx(config['containers'], 'single_cell_pipeline')
    ctx.update(docker_ctx)

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('region'),
        value=normal_bam.keys(),
    )

    workflow.transform(
        name='run_museq',
        ctx=dict(mem=config["memory"]['med'],
                 pool_id=config['pools']['highmem'],
                 **ctx),
        axes=('region',),
        func="single_cell.workflows.mutationseq.tasks.run_museq",
        args=(
            mgd.InputFile("merged_bam", "region", fnames=tumour_bam),
            mgd.InputFile("merged_bai", "region", fnames=tumour_bai),
            mgd.InputFile("normal.split.bam", "region", fnames=normal_bam),
            mgd.InputFile("normal.split.bam.bai", "region", fnames=normal_bai),
            mgd.TempOutputFile("museq.vcf", "region"),
            mgd.TempOutputFile("museq.log", "region"),
            mgd.InputInstance("region"),
            config
        ),
        kwargs={'docker_kwargs': helpers.get_container_ctx(config['containers'], 'mutationseq')}
    )

    workflow.transform(
        name='merge_snvs',
        ctx=dict(mem=config["memory"]['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        func="single_cell.workflows.mutationseq.tasks.concatenate_vcfs",
        args=(
            mgd.TempInputFile("museq.vcf", "region"),
            mgd.TempOutputFile("museq.vcf"),
        ),
    )

    workflow.transform(
        name='finalise_snvs',
        ctx=dict(mem=config["memory"]['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        func="biowrappers.components.io.vcf.tasks.finalise_vcf",
        args=(
            pypeliner.managed.TempInputFile('museq.vcf'),
            pypeliner.managed.OutputFile(snv_vcf),
        ),
        kwargs={'docker_config': helpers.get_container_ctx(config['containers'], 'vcftools')}
    )

    return workflow

