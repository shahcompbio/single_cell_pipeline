'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
import tasks

def create_museq_workflow(normal_bam, normal_bai, tumour_bam, tumour_bai, ref_genome, snv_vcf, snv_csv,
                          config, regions):

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('region'),
        value=regions,
    )


    workflow.transform(name='run_museq',
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

    workflow.transform(name='concatenate_vcfs',
                         ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
                         func=tasks.concatenate_vcfs,
                         args=(
                               mgd.TempInputFile("museq.vcf", "region"),
                               mgd.OutputFile(snv_vcf)
        ),
    )

    workflow.transform(
                         name='parse_museq',
                         ctx={'mem': config["memory"]['low'], 'pool_id': config['pools']['standard'], 'ncpus':1},
                         func=tasks.parse_museq,
                         args=(
                               mgd.InputFile(snv_vcf),
                               mgd.OutputFile(snv_csv),
                               )
                         )

    return workflow

