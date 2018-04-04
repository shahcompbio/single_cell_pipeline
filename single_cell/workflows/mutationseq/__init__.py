'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
import tasks

def create_museq_workflow(tumour_bam, tumour_bai, normal_bam, normal_bai, ref_genome, snv_vcf, snv_csv,
                          config, args, regions):

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('regions'),
        value=regions,
    )


    workflow.transform(name='run_museq',
                         ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['multicore'], 'ncpus':8},
                         axes=('regions',),
                         func=tasks.run_museq,
                         args=(
                               mgd.InputFile("merged_bam", "regions", fnames=tumour_bam),
                               mgd.InputFile("merged_bai", "regions", fnames=tumour_bai),
                               mgd.InputFile(normal_bam),
                               mgd.InputFile(normal_bai),
                               mgd.TempOutputFile("museq.vcf", "regions"),
                               mgd.TempOutputFile("museq.log", "regions"),
                               config,
                               mgd.InputInstance("regions"),
                               regions,
        ),
    )

    workflow.transform(name='concatenate_vcfs',
                         ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['multicore'], 'ncpus':8},
                         func=tasks.concatenate_vcfs,
                         args=(
                               mgd.TempInputFile("museq.vcf", "regions"),
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

