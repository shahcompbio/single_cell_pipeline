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
    regions,
    config):


    normal_split_bam = dict([(ival, normal_split_bam[ival])
                         for ival in regions])
    normal_split_bai = dict([(ival, normal_split_bai[ival])
                         for ival in regions])

    one_split_job = config["one_split_job"]


    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('region'),
        value=regions,
    )

    if one_split_job:
        workflow.transform(
            name='split_normal_bam',
            ctx={'mem': config['memory']['low'], 'pool_id': config['pools']['standard'], 'ncpus':config["max_cores"]},
            func=tasks.split_bam_file_one_job,
            args=(
                mgd.InputFile(normal_bam),
                mgd.InputFile(normal_bai),
                mgd.OutputFile("normal.split.bam", "region", fnames=normal_split_bam, axes_origin=[]),
                mgd.OutputFile("normal.split.bam.bai", "region", fnames=normal_split_bai, axes_origin=[]),
                regions
            ),
            kwargs = {"ncores": config["max_cores"]}
        )

    else:
        workflow.transform(
            name='split_normal_bam',
            ctx={'mem': config['memory']['low']},
            axes=('region',),
            func=tasks.split_bam_file,
            args=(
                mgd.InputFile(normal_bam),
                mgd.InputFile(normal_bai),
                mgd.OutputFile("normal.split.bam", "region", fnames=normal_split_bam),
                mgd.OutputFile("normal.split.bam.bai", "region", fnames=normal_split_bai),
                mgd.InputInstance('region'),
            )
        )

    return workflow
