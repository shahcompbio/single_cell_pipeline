'''
Created on Jul 6, 2017

@author: dgrewal
'''

import pypeliner
import pypeliner.managed as mgd

def create_alignment_workflow(
        bam_filename,
        gc_metrics,
        config,
        ref_genome,
        cell_ids,
):
    baseimage = config['docker']['single_cell_pipeline']

    bam_filename = dict([(cellid, bam_filename[cellid])
                         for cellid in cell_ids])

    workflow = pypeliner.workflow.Workflow(ctx={'docker_image': baseimage})

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=bam_filename.keys(),
    )

    workflow.transform(
        name='bam_collect_gc_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        func="single_cell.utils.picardutils.bam_collect_gc_metrics",
        axes=('cell_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            ref_genome,
            mgd.TempOutputFile('gc_metrics_percell', 'cell_id'),
            mgd.TempOutputFile('gc_metrics_summary_percell', 'cell_id'),
            mgd.TempOutputFile('gc_metrics_pdf_percell', 'cell_id'),
            mgd.TempSpace('gc_tempdir', 'cell_id'),
        ),
        kwargs={
            'docker_image': config['docker']['picard']
        }
    )

    workflow.transform(
        name="collect_gc_metrics",
        func="single_cell.workflows.align.tasks.collect_gc",
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        args=(
            mgd.TempInputFile('gc_metrics_percell', 'cell_id', axes_origin=[]),
            mgd.OutputFile(gc_metrics, extensions=['.yaml']),
            mgd.TempSpace("temp_gc")
        ),
    )

    return workflow
