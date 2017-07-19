'''
Created on Jul 17, 2017

@author: dgrewal
'''


import os
import pypeliner
import pypeliner.managed as mgd
import strelka

import tasks


def create_varcall_workflow(tumour_bam, normal_bam, ref_genome, snv_csv,
                            config, args):
    
    workflow = pypeliner.workflow.Workflow()

    scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
    strelka_parse_script_path = os.path.join(scripts_directory, 'parse_strelka.py')
    museq_parse_script_path = os.path.join(scripts_directory, 'parse_museq.py')
    merge_tables_script = os.path.join(scripts_directory, 'merge.py')


    museq_out_path = os.path.join(args['out_dir'],'pseudo_wgs',
                                  'variant_calling', '{chrom}.mutationseq.vcf')

    museq_out_path_merged = os.path.join(args['out_dir'],'pseudo_wgs',
                                  'variant_calling', 'mutationseq.vcf')

    museq_log_path = os.path.join(args['out_dir'],'pseudo_wgs',
                                  'variant_calling', '{chrom}.mutationseq.log')

    strelka_parsed = os.path.join(args['out_dir'], 'pseudo_wgs', 'variant_calling', 'strelka.csv')
    museq_parsed = os.path.join(args['out_dir'], 'pseudo_wgs', 'variant_calling', 'museq.csv')



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
                             mgd.OutputFile(museq_out_path_merged)
                             )  
                       )

    strelka_indel = os.path.join(args['out_dir'], 'pseudo_wgs', 'variant_calling', 'strelka.indels.vcf.gz')
    strelka_snv = os.path.join(args['out_dir'], 'pseudo_wgs', 'variant_calling', 'strelka.snv.vcf.gz')
    workflow.subworkflow(
            name='strelka',
            func=strelka.create_strelka_workflow,
            args=(
                  mgd.InputFile(args['matched_normal']),
                  mgd.InputFile(tumour_bam),
                  mgd.InputFile(config['ref_genome']),
                  mgd.OutputFile(strelka_indel),
                  mgd.OutputFile(strelka_snv),
            ),
        )

    workflow.commandline(
                         name='parse_strelka',
                         ctx={'mem': config['med_mem']},
                         args=(
                               config['python'],
                               strelka_parse_script_path,
                               '--infile', mgd.InputFile(strelka_snv),
                               '--output', mgd.OutputFile(strelka_parsed),
                               '--tumour_id', 'NA',
                               '--normal_id', 'NA',
                               '--keep_dbsnp','--keep_1000gen',
                               '--remove_duplicates'
                               )
                         )

    workflow.commandline(
                         name='parse_museq',
                         ctx={'mem': config['med_mem']},
                         args=(
                               config['python'],
                               museq_parse_script_path,
                               '--infile', mgd.InputFile(museq_out_path_merged),
                               '--output', mgd.OutputFile(museq_parsed),
                               '--tumour_id', 'NA',
                               '--normal_id', 'NA',
                               '--keep_dbsnp','--keep_1000gen',
                               '--remove_duplicates'
                               )
                         )


    workflow.commandline(
                       name='overlap_var_calls',
                       ctx={'mem': config['med_mem']},
                       args=(
                            config['python'],
                            merge_tables_script,
                            '--merge_type', 'inner',
                            '--nan_value', 'NA',
                            '--input', mgd.InputFile(museq_parsed), mgd.InputFile(strelka_parsed),
                            '--key_cols', 'case_id', 'chromosome', 'start', 'stop', 'ref', 'alt',
                            '--separator', 'tab',
                            '--type', 'merge',
                            '--output', mgd.OutputFile(snv_csv),
                             )
                       )

    return workflow
