'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
import biowrappers.components.io.vcf.tasks
import tasks

def create_museq_workflow(
        normal_bam, normal_bai, tumour_bam, tumour_bai, ref_genome, snv_vcf,
        config):

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('region'),
        value=normal_bam.keys(),
    )

    workflow.transform(
        name='run_museq',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        axes=('region',),
        func=tasks.run_museq,
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
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.concatenate_vcfs,
        args=(
            mgd.TempInputFile("museq.vcf", "region"),
            mgd.TempOutputFile("museq.vcf"),
        ),
    )

    workflow.transform(
        name='finalise_snvs',
        ctx={'pool_id': config['pools']['standard'], 'ncpus':1},
        func=biowrappers.components.io.vcf.tasks.finalise_vcf,
        args=(
            pypeliner.managed.TempInputFile('museq.vcf'),
            pypeliner.managed.OutputFile(snv_vcf),
        ),
    )

    return workflow

