'''
Created on Feb 22, 2018

@author: dgrewal
'''

import os
import pypeliner
import pypeliner.managed as mgd
from workflows import snv_postprocessing
from workflows import germline
from workflows import split_bams
from single_cell.utils import helpers


def germline_calling_workflow(workflow, args):

    config = helpers.load_config(args)

    bam_files, bai_files  = helpers.get_bams(args['input_yaml'])

    sampleids = helpers.get_samples(args['input_yaml'])

    normal_bam_template = args["input_template"]
    normal_bai_template = args["input_template"] + ".bai"

    varcalls_dir = os.path.join(args['out_dir'], 'results',
                                'germline_calling')

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sampleids,
    )
 
    workflow.transform(
        name="get_regions",
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 },
        func=helpers.get_bam_regions,
        ret=pypeliner.managed.TempOutputObj('allregions'),
        args=(
              config["ref_genome"],
              config["split_size"],
              config["chromosomes"],
        )
    )


    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('region'),
        value=mgd.TempInputObj('allregions'),
    )

 
    samtools_germline_vcf = os.path.join(varcalls_dir, 'samtools_germline.vcf.gz')
    samtools_germline_csv = os.path.join(varcalls_dir, 'samtools_germline.csv')
    workflow.subworkflow(
        name='samtools_germline',
        func=germline.create_samtools_germline_workflow,
        args=(
            mgd.InputFile("normal.split.bam", "region", template=normal_bam_template),
            mgd.InputFile("normal.split.bam.bai", "region", template=normal_bai_template),
            config['ref_genome'],
            mgd.OutputFile(samtools_germline_vcf),
            mgd.OutputFile(samtools_germline_csv),
            config,
            mgd.TempInputObj('allregions')
        ),
    )
 
    countdata = os.path.join(varcalls_dir, 'counts.csv')
    olp_calls = os.path.join(varcalls_dir, 'overlapping_calls.csv')
    workflow.subworkflow(
        name='germline_postprocessing',
        func=snv_postprocessing.create_snv_postprocessing_workflow,
        args=(
            mgd.InputFile('bam_markdups', 'sample_id', fnames=bam_files, axes_origin=[]),
            mgd.InputFile('bam_markdups_index', 'sample_id', fnames=bai_files, axes_origin=[]),
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
