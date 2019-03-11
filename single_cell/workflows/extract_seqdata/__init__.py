'''
Created on Apr 13, 2018

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers


def create_extract_seqdata_workflow(
     bam_filename,
     seqdata_filename,
     remixt_config,
     remixt_ref_data_dir,
     config,
):

    ctx = {'mem_retry_increment': 2,
           'docker_image': config['docker']['single_cell_pipeline'],
           'mem': config["memory"]['high']}

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='create_chromosome_seqdata',
        ctx=ctx,
        func="single_cell.workflows.extract_seqdata.tasks.create_chromosome_seqdata",
        args=(
            mgd.OutputFile(seqdata_filename),
            mgd.InputFile(bam_filename, extensions=['.bai']),
            mgd.TempSpace("extract_seqdata_temp"),
            remixt_config,
            remixt_ref_data_dir,
        ),
        kwargs={'chromosomes': config['chromosomes']}
    )

    return workflow
