'''
Created on Nov 21, 2017

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd

import tasks


def create_split_workflow(
    tumour_bam,
    tumour_bai,
    normal_bam,
    normal_bai,
    tumour_split_bam,
    tumour_split_bai,
    normal_split_bam,
    normal_split_bai,
    intervals,
    ref_genome,
    config):


    tumour_split_bam = dict([(ival, tumour_split_bam[ival])
                         for ival in intervals])
    tumour_split_bai = dict([(ival, tumour_split_bai[ival])
                         for ival in intervals])
 
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
            name='split_tumour_bam',
            ctx={'mem': config['memory']['low']},
            func=tasks.split_bam_file_one_job,
            args=(
                mgd.InputFile(tumour_bam),
                mgd.InputFile(tumour_bai),
                ref_genome,
                mgd.OutputFile("tumour.split.bam", "interval", fnames=tumour_split_bam, axes_origin=[]),
                mgd.OutputFile("tumour.split.bam.bai", "interval", fnames=tumour_split_bai, axes_origin=[]),
                intervals
            )
        )

        workflow.transform(
            name='split_normal_bam',
            ctx={'mem': config['memory']['low']},
            func=tasks.split_bam_file_one_job,
            args=(
                mgd.InputFile(normal_bam),
                mgd.InputFile(normal_bai),
                ref_genome,
                mgd.OutputFile("normal.split.bam", "interval", fnames=normal_split_bam, axes_origin=[]),
                mgd.OutputFile("normal.split.bam.bai", "interval", fnames=normal_split_bai, axes_origin=[]),
                intervals
            )
        )

    else:
        workflow.transform(
            name='split_tumour_bam',
            ctx={'mem': config['memory']['low']},
            func=tasks.split_bam_file,
            axes=('interval',),
            args=(
                mgd.InputFile(tumour_bam),
                mgd.InputFile(tumour_bai),
                ref_genome,
                mgd.OutputFile("tumour.split.bam", "interval", fnames=tumour_split_bam),
                mgd.OutputFile("tumour.split.bam.bai", "interval", fnames=tumour_split_bai),
                mgd.InputInstance('interval')
            )
        )

        workflow.transform(
            name='split_normal_bam',
            ctx={'mem': config['memory']['low']},
            axes=('interval',),
            func=tasks.split_bam_file,
            args=(
                mgd.InputFile(normal_bam),
                mgd.InputFile(normal_bai),
                ref_genome,
                mgd.OutputFile("normal.split.bam", "interval", fnames=normal_split_bam),
                mgd.OutputFile("normal.split.bam.bai", "interval", fnames=normal_split_bai),
                mgd.InputInstance('interval')
            )
        )

    return workflow
