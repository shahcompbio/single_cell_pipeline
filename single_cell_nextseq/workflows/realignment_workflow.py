'''
Created on Jul 11, 2017

@author: dgrewal
'''


'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd

import tasks


def create_realignment_workflow(input_bams, output_bams, config, args, sample_ids):

    realn_dir = os.path.join(args['out_dir'], 'realignment')
    targets = os.path.join(realn_dir, 'realignment_targets.intervals')

    output_bams = dict([(sampid, output_bams[sampid]) for sampid in sample_ids])

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )


    if args['realign']:
        workflow.transform(
            name='realignment_targets',
            ctx={'mem': config['high_mem']},
            func=tasks.realigner_target_creator,
            args=(
                mgd.InputFile('bam', 'sample_id', fnames=input_bams),
                mgd.OutputFile(targets),
                mgd.InputFile(config['ref_genome']),
                config
            )
        )
    
        workflow.transform(
            name='gatk_realign',
            ctx={'mem': config['high_mem']},
            func=tasks.gatk_realign,
            args=(
                mgd.InputFile('bam', 'sample_id', fnames=input_bams),
                mgd.OutputFile('bam_realn','sample_id', fnames=output_bams, axes_origin=[]),
                mgd.InputFile(targets),
                mgd.InputFile(config['ref_genome']),
                config,
                mgd.TempSpace('realignment_temp'),
            ),
        )


    else:
        workflow.transform(
                           name='copy_realign',
                           axes=('sample_id',),
                           func = tasks.copy_files,
                           args = (               
                                   mgd.InputFile('bam', 'sample_id', fnames=input_bams),
                                   mgd.OutputFile('bam_copy', 'sample_id', fnames=output_bams),
                                   )
                           )
    
    return workflow
