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
    regions, config, meta_yaml, by_reads=False
):

    ctx = {'mem_retry_increment': 2}
    docker_ctx = helpers.get_container_ctx(config['containers'], 'single_cell_pipeline')
    ctx.update(docker_ctx)

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
            ctx=dict(mem=config['memory']['low'],
                     pool_id=config['pools']['multicore'],
                     ncpus=config['max_cores'],
                     **ctx),
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
                helpers.get_container_ctx(config['containers'], 'samtools')
            ),
        )

    elif one_split_job:
        workflow.transform(
            name='split_normal_bam',
            ctx=dict(mem=config['memory']['low'],
                     pool_id=config['pools']['multicore'],
                     ncpus=config['max_cores'],
                     **ctx),
            func="single_cell.workflows.split_bams.tasks.split_bam_file_one_job",
            args=(
                mgd.InputFile(normal_bam, extensions=['.bai']),
                mgd.OutputFile(
                    "normal.split.bam", "region",
                    fnames=normal_split_bam, axes_origin=[],
                    extensions=['.bai'],
                ),
                regions,
                helpers.get_container_ctx(config['containers'], 'samtools')
            ),
            kwargs={"ncores": config["max_cores"]}
        )

    else:
        workflow.transform(
            name='split_normal_bam',
            ctx=dict(mem=config['memory']['low'],
                     pool_id=config['pools']['standard'],
                     ncpus=1,
                     **ctx),
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
                helpers.get_container_ctx(config['containers'], 'samtools')
            )
        )


    inputs = helpers.format_file_yaml(normal_bam)
    outputs = {k: helpers.format_file_yaml(v) for k,v in normal_split_bam.iteritems()}

    metadata = {
        'split_bam': {
            'name': 'split_bam',
            'ref_genome': config["ref_genome"],
            'version': single_cell.__version__,
            'containers': config['containers'],
            'output_datasets': outputs,
            'input_datasets': inputs,
            'split_by_reads': by_reads,
            'one_job_for_splitting': one_split_job
        }
    }

    workflow.transform(
        name='generate_meta_yaml',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        func="single_cell.utils.helpers.write_to_yaml",
        args=(
            mgd.OutputFile(meta_yaml),
            metadata
        )
    )



    return workflow
