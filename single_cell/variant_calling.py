'''
Created on Feb 22, 2018

@author: dgrewal
'''

import os
import pypeliner
import pypeliner.managed as mgd
from workflows import snv_postprocessing
from workflows import mutationseq 
from workflows import strelka
from single_cell.utils import helpers

def variant_calling_workflow(workflow, args):

    config = helpers.load_config(args)

    pseudo_wgs_bam = os.path.join(args['out_dir'], 'pseudo_wgs',
                                  'merged.sorted.markdups.bam')
    pseudo_wgs_bai = os.path.join(args['out_dir'], 'pseudo_wgs',
                                  'merged.sorted.markdups.bam.bai')

    bam_files, bai_files  = helpers.get_bams(args['bams_file'])

    sampleids = helpers.get_samples(args['bams_file'])

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sampleids,
    )

    workflow.transform(
        name='generate_intervals',
        func=helpers.generate_intervals,
        ret=pypeliner.managed.TempOutputObj('intervals'),
        args=(
            config["ref_genome"],
        )
    )
    
    varcalls_dir = os.path.join(args['out_dir'], 'pseudo_wgs',
                                'variant_calling')
    museq_vcf = os.path.join(varcalls_dir, 'museq_snv.vcf')
    museq_csv = os.path.join(varcalls_dir, 'museq_snv.csv')
    workflow.subworkflow(
        name='museq',
        func=mutationseq.create_museq_workflow,
        args=(
            mgd.InputFile(pseudo_wgs_bam),
            mgd.InputFile(pseudo_wgs_bai),
            mgd.InputFile(args['matched_normal']),
            mgd.InputFile(args['matched_normal'] + ".bai"),
            config['ref_genome'],
            mgd.OutputFile(museq_vcf),
            mgd.OutputFile(museq_csv),
            config,
            args,
            mgd.TempInputObj('intervals')
        ),
    )
    
    strelka_snv_vcf = os.path.join(varcalls_dir, 'strelka_snv.vcf')
    strelka_indel_vcf = os.path.join(varcalls_dir, 'strelka_indel.vcf')
    strelka_snv_csv = os.path.join(varcalls_dir, 'strelka_snv.csv')
    strelka_indel_csv = os.path.join(varcalls_dir, 'strelka_indel.csv')
    workflow.subworkflow(
        name='strelka',
        func=strelka.create_strelka_workflow,
        args=(
            mgd.InputFile(args['matched_normal']),
            mgd.InputFile(args['matched_normal'] + ".bai"),
            mgd.InputFile(pseudo_wgs_bam),
            mgd.InputFile(pseudo_wgs_bai),
            config['ref_genome'],
            mgd.OutputFile(strelka_indel_vcf),
            mgd.OutputFile(strelka_snv_vcf),
            mgd.OutputFile(strelka_indel_csv),
            mgd.OutputFile(strelka_snv_csv),
            config,
            mgd.TempInputObj('intervals')
        ),
    )
     
    countdata = os.path.join(args['out_dir'], 'pseudo_wgs',
                             'counts', 'counts.csv')
    workflow.subworkflow(
        name='postprocessing',
        func=snv_postprocessing.create_snv_postprocessing_workflow,
        args=(
            mgd.InputFile('bam_markdups', 'sample_id', fnames=bam_files),
            mgd.InputFile('bam_markdups_index', 'sample_id', fnames=bai_files),
            mgd.InputFile(museq_csv),
            mgd.InputFile(strelka_snv_csv),
            mgd.OutputFile(countdata),
            sampleids,
            config,
            args['out_dir'],
        )
    )

    return workflow
