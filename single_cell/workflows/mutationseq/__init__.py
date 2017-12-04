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

    tumour_bam = dict([(ival, tumour_bam[ival])
                         for ival in intervals])
  
    normal_bam = dict([(ival, normal_bam[ival])
                         for ival in intervals])

    workflow.setobj(
        obj=mgd.OutputChunks('interval'),
        value=intervals,
    )

    workflow.transform(name='run_museq',
                         ctx={'mem': config['med_mem']},
                         axes = ('interval',),
                         func=tasks.run_museq,
                         args=(
                               mgd.InputFile("tumour.split.bam", "interval", fnames=tumour_bam),
                               mgd.InputFile("tumour.split.bam.bai", "interval", fnames=tumour_bai),
                               mgd.InputFile("normal.split.bam", "interval", fnames=normal_bam),
                               mgd.InputFile("normal.split.bam.bai", "interval", fnames=normal_bai),
                               config['ref_genome'],
                               config['mutationseq'],
                               mgd.TempOutputFile("museq.vcf", 'interval'),
                               mgd.TempOutputFile("museq.log", 'interval'),
                               mgd.InputInstance('interval'),
                               config
                               )
                         )

    workflow.transform(
                       name='merge_vcf',
                       ctx={'mem': config['low_mem']},
                       func=tasks.concatenate_vcf,
                       args=(mgd.TempInputFile("museq.vcf", 'interval'),
                             mgd.OutputFile(snv_vcf)
                             )  
                       )


    workflow.transform(
                         name='parse_museq',
                         ctx={'mem': config['low_mem']},
                         func=tasks.parse_museq,
                         args=(
                               mgd.InputFile(snv_vcf),
                               mgd.OutputFile(snv_csv),
                               )
                         )

    return workflow

