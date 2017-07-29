'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd

import tasks


def create_realignment_workflow(input_bams, output_bams, config, out_dir, realign, sample_ids):

    realn_dir = os.path.join(out_dir, 'realignment')
    targets = os.path.join(realn_dir, 'realignment_targets.intervals')

    output_bams = dict([(sampid, output_bams[sampid]) for sampid in sample_ids])

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )

    if realign:

        chromosomes = map(str, range(1,22)) + ['X', 'Y']
    
        workflow.setobj(
            obj=mgd.OutputChunks('chrom'),
            value=chromosomes,
        )

        workflow.transform(
            name='realignment_targets',
            axes=('chrom',),
            ctx={'mem': config['high_mem']},
            func=tasks.realign,
            args=(
                mgd.InputFile('bam', 'sample_id', fnames=input_bams),
                mgd.TempOutputFile('realigned.bam', 'chrom', 'sample_id'),
                mgd.TempSpace('realignment_temp', 'chrom', cleanup='before'),
                config,
                mgd.InputInstance('chrom')
            )
        )


        workflow.transform(
            name='merge_realignment',
            ctx={'mem': config['high_mem']},
            axes=('sample_id',),
            func=tasks.merge_realignment,
            args=(
                  mgd.TempInputFile('realigned.bam', 'chrom', 'sample_id'),
                  mgd.OutputFile('bam_realn','sample_id', fnames=output_bams),
                  config,
                  mgd.InputInstance('sample_id')
                  )
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
