'''
Created on Nov 21, 2017

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd

import tasks


def create_split_workflow(
    normal_bam,
    normal_bai,
    normal_split_bam,
    normal_split_bai,
    intervals,
    config):


    normal_split_bam = dict([(ival, normal_split_bam[ival])
                         for ival in intervals])
    normal_split_bai = dict([(ival, normal_split_bai[ival])
                         for ival in intervals])

    one_split_job = config["one_split_job"]


    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('interval'),
        value=intervals,
    )

    if one_split_job:
        workflow.transform(
            name='split_normal_bam',
            ctx={'mem': config['memory']['low']},
            func=tasks.split_bam_file_one_job,
            args=(
                mgd.InputFile(normal_bam),
                mgd.InputFile(normal_bai),
                mgd.OutputFile("normal.split.bam", "interval", fnames=normal_split_bam, axes_origin=[]),
                mgd.OutputFile("normal.split.bam.bai", "interval", fnames=normal_split_bai, axes_origin=[]),
                intervals
            )
        )

    else:
        workflow.transform(
            name='split_normal_bam',
            ctx={'mem': config['memory']['low']},
            axes=('interval',),
            func=tasks.split_bam_file,
            args=(
                mgd.InputFile(normal_bam),
                mgd.InputFile(normal_bai),
                mgd.OutputFile("normal.split.bam", "interval", fnames=normal_split_bam),
                mgd.OutputFile("normal.split.bam.bai", "interval", fnames=normal_split_bai),
                mgd.InputInstance('interval')
            )
        )

    return workflow
