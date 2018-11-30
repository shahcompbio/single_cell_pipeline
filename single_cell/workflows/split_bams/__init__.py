'''
Created on Nov 21, 2017

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers
import single_cell

def create_split_workflow(
    normal_bam, normal_bai, normal_split_bam, normal_split_bai,
    regions, config, by_reads=False
):

    baseimage = config['docker']['single_cell_pipeline']
    samtoolsimage = config['docker']['samtools']

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

    # split by reads always runs no a single node
    if by_reads:
        workflow.transform(
            name='split_normal_bam',
            ctx={'mem': config['memory']['low'], 'ncpus': config['max_cores'], 'docker_image': baseimage},
            func="single_cell.workflows.split_bams.tasks.split_bam_file_by_reads",
            args=(
                mgd.InputFile(normal_bam),
                mgd.InputFile(normal_bai),
                mgd.OutputFile(
                    "normal.split.bam", "region",
                    fnames=normal_split_bam, axes_origin=[]
                ),
                mgd.OutputFile(
                    "normal.split.bam.bai", "region",
                    fnames=normal_split_bai, axes_origin=[]
                ),
                mgd.TempSpace("bam_split_by_reads"),
                regions,
                samtoolsimage
            ),
        )

    elif one_split_job:
        workflow.transform(
            name='split_normal_bam',
            ctx={'mem': config['memory']['low'], 'ncpus': config['max_cores'], 'docker_image': baseimage},
            func="single_cell.workflows.split_bams.tasks.split_bam_file_one_job",
            args=(
                mgd.InputFile(normal_bam, extensions=['.bai']),
                mgd.OutputFile(
                    "normal.split.bam", "region",
                    fnames=normal_split_bam, axes_origin=[],
                    extensions=['.bai'],
                ),
                regions,
                samtoolsimage,
                mgd.TempSpace("one_job_split_tempdir")
            ),
            kwargs={"ncores": config["max_cores"]}
        )

    else:
        workflow.transform(
            name='split_normal_bam',
            ctx={'mem': config['memory']['low'], 'ncpus': config['max_cores'], 'docker_image': baseimage},
            axes=('region',),
            func="single_cell.workflows.split_bams.tasks.split_bam_file",
            args=(
                mgd.InputFile(normal_bam),
                mgd.InputFile(normal_bai),
                mgd.OutputFile(
                    "normal.split.bam", "region", fnames=normal_split_bam
                ),
                mgd.OutputFile(
                    "normal.split.bam.bai", "region", fnames=normal_split_bai
                ),
                mgd.InputInstance('region'),
                samtoolsimage
            )
        )


    return workflow
