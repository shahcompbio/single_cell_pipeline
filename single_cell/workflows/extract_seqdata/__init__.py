'''
Created on Apr 13, 2018

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd

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

    singlecellimage = config['docker']['images']['single_cell_pipeline']
    ctx = {
              'mem_retry_increment': 2,
              'ncpus': 1,
              'image': singlecellimage['image'],
              'dockerize': config['docker']['dockerize'],
              'mounts': config['docker']['mounts'],
              'username': singlecellimage['username'],
              'password': singlecellimage['password'],
              'server': singlecellimage['server'],
          }

    workflow = pypeliner.workflow.Workflow()


    workflow.transform(
        name="get_chromosomes",
        ctx=dict(mem=2, pool_id=config['pools']['standard'], **ctx),
        func="remixt.config.get_chromosomes",
        ret=mgd.TempOutputObj('chromosomes'),
        args=(
              remixt_config,
              ref_data_dir
        )
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

