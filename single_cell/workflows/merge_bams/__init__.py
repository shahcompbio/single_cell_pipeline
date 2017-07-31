'''
Created on Jul 11, 2017

@author: dgrewal
'''


'''
Created on Jul 6, 2017

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd

import tasks


def create_merge_workflow(
    bam,
    bam_filename,
    bam_index_filename,
    config,
    lanes):
 
    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('lane'),
        value=lanes,
    )


    workflow.transform(
        name='merge_bams',
        ctx={'mem': config['high_mem']},
        func=tasks.merge_bams,
        args=(
            mgd.InputFile('bam', 'lane', fnames=bam),
            mgd.OutputFile(bam_filename),
        ),
    )

    workflow.commandline(
        name='bam_index',
        ctx={'mem': config['low_mem']},
        args=(
            'samtools', 'index',
            mgd.InputFile(bam_filename),
            mgd.OutputFile(bam_index_filename),
        ),
    )



    return workflow
