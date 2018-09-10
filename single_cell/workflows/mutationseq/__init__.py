'''
Created on Jul 24, 2017

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd

def create_museq_workflow(
        normal_bam, normal_bai, tumour_bam, tumour_bai, ref_genome, snv_vcf,
        config):

    singlecellimage = config['docker']['images']['single_cell_pipeline']

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('region'),
        value=normal_bam.keys(),
    )

    workflow.transform(
        name='run_museq',
        ctx={
            'mem': config["memory"]['med'], 'num_retry': 3,
            'mem_retry_increment': 2,
            'pool_id': config['pools']['highmem'], 'ncpus': 1,
            'image': singlecellimage['image'],
            'dockerize': config['docker']['dockerize'],
            'mounts': config['docker']['mounts'],
            'username': singlecellimage['username'],
            'password': singlecellimage['password'],
            'server': singlecellimage['server'],
        },
        axes=('region',),
        func="single_cell.workflows.mutationseq.tasks.run_museq",
        args=(
            mgd.InputFile("merged_bam", "region", fnames=tumour_bam),
            mgd.InputFile("merged_bai", "region", fnames=tumour_bai),
            mgd.InputFile("normal.split.bam", "region", fnames=normal_bam),
            mgd.InputFile("normal.split.bam.bai", "region", fnames=normal_bai),
            mgd.TempOutputFile("museq.vcf", "region"),
            mgd.TempOutputFile("museq.log", "region"),
            config,
            mgd.InputInstance("region"),
        ),
    )

    workflow.transform(
        name='merge_snvs',
        ctx={
            'mem': config["memory"]['med'], 'num_retry': 3,
            'mem_retry_increment': 2,
            'pool_id': config['pools']['standard'], 'ncpus': 1,
            'image': singlecellimage['image'],
            'dockerize': config['docker']['dockerize'],
            'mounts': config['docker']['mounts'],
            'username': singlecellimage['username'],
            'password': singlecellimage['password'],
            'server': singlecellimage['server'],
        },
        func="single_cell.workflows.mutationseq.tasks.concatenate_vcfs",
        args=(
            mgd.TempInputFile("museq.vcf", "region"),
            mgd.TempOutputFile("museq.vcf"),
        ),
    )

    workflow.transform(
        name='finalise_snvs',
        ctx={
            'mem': config["memory"]['med'], 'num_retry': 3,
            'mem_retry_increment': 2,
            'pool_id': config['pools']['standard'], 'ncpus': 1,
            'image': singlecellimage['image'],
            'dockerize': config['docker']['dockerize'],
            'mounts': config['docker']['mounts'],
            'username': singlecellimage['username'],
            'password': singlecellimage['password'],
            'server': singlecellimage['server'],
        },
        func="biowrappers.components.io.vcf.tasks.finalise_vcf",
        args=(
            pypeliner.managed.TempInputFile('museq.vcf'),
            pypeliner.managed.OutputFile(snv_vcf),
        ),
        kwargs={'docker_config': config['docker']}
    )

    return workflow

