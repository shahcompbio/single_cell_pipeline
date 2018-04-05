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


def germline_calling_workflow(workflow, args):

    config = helpers.load_config(args)

    bam_files, bai_files  = helpers.get_bams(args['input_yaml'])

    sampleids = helpers.get_samples(args['input_yaml'])

    varcalls_dir = os.path.join(args['out_dir'], 'results',
                                'germline_calling')

    wgs_bam_dir = args["merged_wgs"]
    wgs_bam_template = os.path.join(wgs_bam_dir, "{regions}_merged.bam")
    wgs_bai_template = os.path.join(wgs_bam_dir, "{regions}_merged.bam.bai")

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sampleids,
    )

    workflow.transform(
        name="get_regions",
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 },
        func=helpers.get_bam_regions,
        ret=pypeliner.managed.TempOutputObj('regions'),
        args=(
              config["ref_genome"],
              config["split_size"],
              config["chromosomes"],
        )
    )

    workflow.subworkflow(
        name="split_normal",
        func=split_bams.create_split_workflow,
        args = (
            mgd.InputFile(args['matched_normal']),
            mgd.InputFile(args['matched_normal'] + ".bai"),
            mgd.TempOutputFile("normal.split.bam", "regions", axes_origin=[]),
            mgd.TempOutputFile("normal.split.bam.bai", "regions", axes_origin=[]),
            pypeliner.managed.TempInputObj('regions'),
            config,
        )
    )

    samtools_germline_vcf = os.path.join(varcalls_dir, 'samtools_germline.vcf')
    samtools_germline_csv = os.path.join(varcalls_dir, 'samtools_germline.csv')
    workflow.subworkflow(
        name='samtools_germline',
        func=germline.create_samtools_germline_workflow,
        args=(
            mgd.TempInputFile("normal.split.bam", "regions"),
            mgd.TempInputFile("normal.split.bam.bai", "regions"),
            config['ref_genome'],
            mgd.OutputFile(samtools_germline_vcf),
            mgd.OutputFile(samtools_germline_csv),
            config,
            mgd.TempInputObj('regions')
        ),
    )

    countdata = os.path.join(varcalls_dir, 'counts.csv')
    olp_calls = os.path.join(varcalls_dir, 'overlapping_calls.csv')
    workflow.subworkflow(
        name='postprocessing',
        func=snv_postprocessing.create_snv_postprocessing_workflow,
        args=(
            mgd.InputFile('bam_markdups', 'sample_id', fnames=bam_files),
            mgd.InputFile('bam_markdups_index', 'sample_id', fnames=bai_files),
            [
                mgd.InputFile(samtools_germline_csv),
            ],
            mgd.OutputFile(countdata),
            mgd.OutputFile(olp_calls),
            sampleids,
            config,
        )
    )

    return workflow
