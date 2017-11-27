'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
import tasks

def create_museq_workflow(tumour_bam, normal_bam, ref_genome, snv_vcf, snv_csv,
                            config, args):
    
    workflow = pypeliner.workflow.Workflow()

    museq_out_path = os.path.join(args['out_dir'],'pseudo_wgs',
                                  'variant_calling', '{chrom}.mutationseq.vcf')

    museq_log_path = os.path.join(args['out_dir'],'pseudo_wgs',
                                  'variant_calling', '{chrom}.mutationseq.log')


    chromosomes = config["chromosomes"]

    workflow.setobj(
        obj=mgd.OutputChunks('chrom'),
        value=chromosomes,
    )

    workflow.transform(name='run_museq',
                         axes=('chrom',),
                         ctx={'mem': config['med_mem']},
                         func=tasks.run_museq,
                         args=(
                               mgd.InputFile(tumour_bam),
                               mgd.InputFile(normal_bam),
                               config['ref_genome'],
                               config['mutationseq'],
                               mgd.OutputFile(museq_out_path, 'chrom'),
                               mgd.OutputFile(museq_log_path, 'chrom'),
                               mgd.InputInstance('chrom'),
                               config
                               )
                         )

    workflow.transform(
                       name='merge_vcf',
                       ctx={'mem': config['low_mem']},
                       func=tasks.concatenate_vcf,
                       args=(mgd.InputFile(museq_out_path, 'chrom'),
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

