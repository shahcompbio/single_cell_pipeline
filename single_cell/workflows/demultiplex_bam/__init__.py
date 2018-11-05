'''
Created on November 01, 2018

@author: pashaa
'''
import os
import pypeliner
from single_cell.utils import helpers


def create_demultiplex_bam_workflow(
        managed_input_file,
        config,
        out_dir,
        barcode_csv):

    ctx = {'mem_retry_increment': 2, 'ncpus': 1}
    docker_ctx = helpers.get_container_ctx(config['containers'], 'single_cell_pipeline')
    ctx.update(docker_ctx)

    workflow = pypeliner.workflow.Workflow()

    # demultiplex
    workflow.transform(
        name='demultiplex',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        func="single_cell.workflows.demultiplex_bam.tasks.demultiplex",
        args=(managed_input_file,
              pypeliner.managed.TempOutputFile('demux', 'split'),
              barcode_csv)
    )

    # collate demultiplexed bams
    workflow.transform(
        name='collate',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        axes=('split',),
        func="single_cell.workflows.demultiplex_bam.tasks.collate",
        args=(pypeliner.managed.TempInputFile('demux', 'split'),
              pypeliner.managed.TempOutputFile('collated', 'split'))
    )

    # convert collated bams to fastq
    outfile_r1 = os.path.join(out_dir, '{split}_r1.fastq.gz')
    outfile_r2 = os.path.join(out_dir, '{split}_r2.fastq.gz')

    workflow.transform(
        name='fastq',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        axes=('split',),
        func="single_cell.workflows.demultiplex_bam.tasks.fastq",
        args=(pypeliner.managed.TempInputFile('collated', 'split'),
              pypeliner.managed.OutputFile('fq1', 'split', template=outfile_r1),
              pypeliner.managed.OutputFile('fq2', 'split', template=outfile_r2)))

    return workflow
