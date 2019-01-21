'''
Created on Jul 24, 2017

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers


def create_museq_workflow(
        normal_bam, tumour_bam, ref_genome, snv_vcf,
        config):


    museq_docker = {'docker_image': config['docker']['mutationseq']}
    vcftools_docker = {'docker_image': config['docker']['vcftools']}

    ctx = {'mem_retry_increment': 2, 'ncpus': 1, 'num_retry': 3,
           'docker_image': config['docker']['single_cell_pipeline']}

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('region'),
        value=normal_bam.keys(),
    )

    workflow.transform(
        name='run_museq',
        ctx=dict(mem=config["memory"]['med'],
                 **ctx),
        axes=('region',),
        func='single_cell.workflows.mutationseq.tasks.run_museq',
        args=(
            mgd.InputFile('merged_bam', 'region', fnames=tumour_bam, extensions=['.bai']),
            mgd.InputFile('normal.split.bam', 'region', fnames=normal_bam, extensions=['.bai']),
            mgd.TempOutputFile('museq.vcf', 'region'),
            mgd.TempOutputFile('museq.log', 'region'),
            mgd.InputInstance('region'),
            config,
        ),
        kwargs={'docker_kwargs': museq_docker}
    )

    workflow.transform(
        name='finalise_region_vcfs',
        axes=('region',),
        ctx=dict(mem=config["memory"]['med'], **ctx),
        func='biowrappers.components.io.vcf.tasks.finalise_vcf',
        args=(
            mgd.TempInputFile('museq.vcf', 'region'),
            mgd.TempOutputFile('museq.vcf.gz', 'region', extensions=['.tbi', '.csi']),
        ),
        kwargs={'docker_config': vcftools_docker}
    )

    workflow.transform(
        name='merge_snvs',
        ctx=dict(mem=config["memory"]['med'],
                 **ctx),
        func='biowrappers.components.io.vcf.tasks.concatenate_vcf',
        args=(
            mgd.TempInputFile('museq.vcf.gz', 'region', extensions=['.tbi', '.csi']),
            mgd.TempOutputFile('museq.vcf.gz', extensions=['.tbi', '.csi']),
        ),
        kwargs={
            'allow_overlap': True,
            'docker_config': vcftools_docker
        },
    )

    workflow.transform(
        name='finalise_vcf',
        func='biowrappers.components.io.vcf.tasks.finalise_vcf',
        ctx=dict(mem=config["memory"]['med'],
                 **ctx),
        args=(
            mgd.TempInputFile('museq.vcf.gz', extensions=['.tbi', '.csi']),
            mgd.OutputFile(snv_vcf, extensions=['.tbi', '.csi']),
        ),
        kwargs={'docker_config': vcftools_docker}
    )

    return workflow

