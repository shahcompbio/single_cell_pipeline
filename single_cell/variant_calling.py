'''
Created on Feb 22, 2018

@author: dgrewal
'''

import os
import pypeliner
import pypeliner.managed as mgd
from workflows import snv_postprocessing
from workflows import mutationseq 
from workflows import split_bams
from workflows import strelka
from single_cell.utils import helpers

def variant_calling_workflow(workflow, args):

    config = helpers.load_config(args)


    bam_files, bai_files  = helpers.get_bams(args['input_yaml'])

    cellids = helpers.get_samples(args['input_yaml'])

    varcalls_dir = os.path.join(args['out_dir'], 'results',
                                'variant_calling')


    wgs_bam_template = args["merged_wgs_template"]
    wgs_bai_template = args["merged_wgs_template"] + ".bai"

    normal_bam_template = args["normal_split_template"]
    normal_bai_template = args["normal_split_template"] + ".bai"


    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cellids,
    )

    workflow.transform(
        name="get_regions",
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 },
        func=helpers.get_bam_regions,
        ret=pypeliner.managed.TempOutputObj('region'),
        args=(
              config["ref_genome"],
              config["split_size"],
              config["chromosomes"],
        )
    )

    museq_vcf = os.path.join(varcalls_dir, 'museq_snv.vcf')
    museq_csv = os.path.join(varcalls_dir, 'museq_snv.csv')
    workflow.subworkflow(
        name='museq',
        func=mutationseq.create_museq_workflow,
        args=(
            mgd.InputFile("normal.split.bam", "region", template=normal_bam_template, axes_origin=[]),
            mgd.InputFile("normal.split.bam.bai", "region", template=normal_bai_template, axes_origin=[]),
            mgd.InputFile("merged_bam", "region", template=wgs_bam_template, axes_origin=[]),
            mgd.InputFile("merged_bam", "region", template=wgs_bai_template, axes_origin=[]),
            config['ref_genome'],
            mgd.OutputFile(museq_vcf),
            mgd.OutputFile(museq_csv),
            config,
            mgd.TempInputObj('region')
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
            mgd.InputFile("normal.split.bam", "region", template=normal_bam_template, axes_origin=[]),
            mgd.InputFile("normal.split.bam.bai", "region", template=normal_bai_template, axes_origin=[]),
            mgd.InputFile("merged_bam", "region", template=wgs_bam_template, axes_origin=[]),
            mgd.InputFile("merged_bam", "region", template=wgs_bai_template, axes_origin=[]),
            config['ref_genome'],
            mgd.OutputFile(strelka_indel_vcf),
            mgd.OutputFile(strelka_snv_vcf),
            mgd.OutputFile(strelka_indel_csv),
            mgd.OutputFile(strelka_snv_csv),
            config,
            mgd.TempInputObj('region')
        ),
    )

    countdata = os.path.join(varcalls_dir, 'counts.csv')
    olp_calls = os.path.join(varcalls_dir, 'overlapping_calls.csv')
    workflow.subworkflow(
        name='postprocessing',
        func=snv_postprocessing.create_snv_postprocessing_workflow,
        args=(
            mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_files),
            mgd.InputFile('bam_markdups_index', 'cell_id', fnames=bai_files),
            [
                mgd.InputFile(museq_csv),
                mgd.InputFile(strelka_snv_csv),
            ],
            mgd.OutputFile(countdata),
            mgd.OutputFile(olp_calls),
            cellids,
            config,
        )
    )

    return workflow
