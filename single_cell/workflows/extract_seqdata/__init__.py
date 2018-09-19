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
     multiprocess=False,
):
    ctx = {'mem_retry_increment': 2, 'ncpus': 1}
    docker_ctx = helpers.get_container_ctx(config['containers'], 'single_cell_pipeline')
    ctx.update(docker_ctx)

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='create_chromosome_seqdata',
        ctx=dict(mem=config["memory"]['high'], pool_id=config['pools']['highmem'], **ctx),
        func="single_cell.workflows.extract_seqdata.tasks.create_chromosome_seqdata",
        args=(
            mgd.TempOutputFile('seqdata', 'chromosome'),
            mgd.InputFile(bam_filename),
            remixt_config,
            remixt_ref_data_dir,
        ),
        kwargs={'multiprocess':multiprocess,
                'ncores':config['max_cores']}
    )

    workflow.transform(
        name='merge_seqdata',
        ctx=dict(mem=config["memory"]['high'], pool_id=config['pools']['highmem'], **ctx),
        func="remixt.seqdataio.merge_seqdata",
        args=(
            mgd.OutputFile(seqdata_filename),
            mgd.TempInputFile('seqdata', 'chromosome'),
        ),
    )

    return workflow
