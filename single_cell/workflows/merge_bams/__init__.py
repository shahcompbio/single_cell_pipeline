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
    sample_lanes):


    lanes = list(set([v[1] for v in sample_lanes]))
    samples = list(set([v[0] for v in sample_lanes]))

    bam_filename = dict([(sampid, bam_filename[sampid])for sampid in samples])
    bam_index_filename = dict([(sampid, bam_index_filename[sampid])for sampid in samples])

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('lane'),
        value=lanes,
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'lane'),
        value=sample_lanes,
    )

    workflow.transform(
        name='merge_bams',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard']},
        func=tasks.merge_bams,
        axes=('sample_id',),
        args=(
            mgd.InputFile('bam', 'sample_id', 'lane', fnames=bam),
            mgd.OutputFile('merged_lanes.bam', 'sample_id', fnames=bam_filename),
            mgd.OutputFile('merged_lanes.bam.bai', 'sample_id', fnames=bam_index_filename),
        ),
    )

    return workflow
