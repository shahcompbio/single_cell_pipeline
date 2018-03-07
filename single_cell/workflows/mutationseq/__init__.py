'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
import tasks

def create_museq_workflow(tumour_bam, tumour_bai, normal_bam, normal_bai, ref_genome, snv_vcf, snv_csv,
                          config, args, intervals):
    
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(name='run_museq',
                         ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['multicore'], 'ncpus':8},
                         func=tasks.run_museq,
                         args=(
                               mgd.InputFile(tumour_bam),
                               mgd.InputFile(tumour_bai),
                               mgd.InputFile(normal_bam),
                               mgd.InputFile(normal_bai),
                               mgd.OutputFile(snv_vcf),
                               mgd.TempOutputFile("museq.log"),
                               mgd.TempSpace("temp"),
                               config,
                               intervals,
        ),
        kwargs={"ncores": config["max_cores"]},
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

