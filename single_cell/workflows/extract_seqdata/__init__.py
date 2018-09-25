'''
Created on Apr 13, 2018

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers


def create_extract_seqdata_workflow(
     bam_filename,
     bai_filename,
     seqdata_filename,
     config,
     remixt_config,
     ref_data_dir,
     snp_positions_filename,
     bam_max_fragment_length,
     bam_max_soft_clipped,
     bam_check_proper_pair,
     multiprocess=False
):


    ctx = {'mem_retry_increment': 2, 'ncpus': 1}
    docker_ctx = helpers.get_container_ctx(config['containers'], 'single_cell_pipeline')
    ctx.update(docker_ctx)

    workflow = pypeliner.workflow.Workflow()


    workflow.transform(
        name="get_chromosomes",
        ctx=dict(mem=2, pool_id=config['pools']['standard'], **ctx),
        func="single_cell.workflows.extract_seqdata.tasks.get_chromosomes",
        ret=mgd.TempOutputObj('chromosomes'),
        args=(
              remixt_config,
              ref_data_dir,
        ),
        kwargs={'chromosomes':config['chromosomes']},
    )

    workflow.setobj(
        obj=mgd.OutputChunks('chromosome'),
        value=mgd.TempInputObj('chromosomes')
    )

    workflow.transform(
        name='create_chromosome_seqdata',
        ctx=dict(mem=config["memory"]['high'], pool_id=config['pools']['highmem'], **ctx),
        func="single_cell.workflows.extract_seqdata.tasks.create_chromosome_seqdata",
        args=(
            mgd.TempOutputFile('seqdata', 'chromosome', axes_origin=[]),
            mgd.InputFile(bam_filename),
            mgd.InputFile(bai_filename),
            snp_positions_filename,
            mgd.TempInputObj('chromosomes'),
            bam_max_fragment_length,
            bam_max_soft_clipped,
            bam_check_proper_pair,
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

