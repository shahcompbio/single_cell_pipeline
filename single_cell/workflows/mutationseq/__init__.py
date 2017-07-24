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

    scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')

    museq_parse_script_path = os.path.join(scripts_directory, 'parse_museq.py')


    museq_out_path = os.path.join(args['out_dir'],'pseudo_wgs',
                                  'variant_calling', '{chrom}.mutationseq.vcf')

    museq_out_path_merged = os.path.join(args['out_dir'],'pseudo_wgs',
                                  'variant_calling', 'mutationseq.vcf')

    museq_log_path = os.path.join(args['out_dir'],'pseudo_wgs',
                                  'variant_calling', '{chrom}.mutationseq.log')


    chromosomes = map(str, range(1,22)) + ['X', 'Y']

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
                               mgd.InputFile(config['ref_genome']),
                               config['mutationseq'],
                               mgd.OutputFile(museq_out_path, 'chrom'),
                               mgd.OutputFile(museq_log_path, 'chrom'),
                               mgd.InputInstance('chrom'),
                               config
                               )
                         )

    workflow.transform(
                       name='merge_vcf',
                       ctx={'mem': config['med_mem']},
                       func=tasks.concatenate_vcf,
                       args=(mgd.InputFile(museq_out_path, 'chrom'),
                             mgd.OutputFile(snv_vcf)
                             )  
                       )


    workflow.commandline(
                         name='parse_museq',
                         ctx={'mem': config['med_mem']},
                         args=(
                               config['python'],
                               museq_parse_script_path,
                               '--infile', mgd.InputFile(museq_out_path_merged),
                               '--output', mgd.OutputFile(snv_csv),
                               '--tumour_id', 'NA',
                               '--normal_id', 'NA',
                               '--keep_dbsnp','--keep_1000gen',
                               '--remove_duplicates'
                               )
                         )

    return workflow

