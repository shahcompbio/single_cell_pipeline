'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
import single_cell
from single_cell.utils import helpers

def create_alignment_workflow(
        fastq_1_filename,
        fastq_2_filename,
        bam_filename,
        alignment_metrics,
        gc_metrics,
        plot_metrics,
        ref_genome,
        config,
        args,
        triminfo,
        centerinfo,
        sample_info,
        cell_ids,
):

    disable_biobloom = args['disable_biobloom']

    baseimage = config['docker']['single_cell_pipeline']

    out_dir = args['out_dir']

    merge_metrics = os.path.join(out_dir, 'metrics')

    lane_metrics = os.path.join(args['out_dir'], 'metrics_per_lane', '{lane}')

    bam_filename = dict([(cellid, bam_filename[cellid])
                         for cellid in cell_ids])

    chromosomes = config["chromosomes"]

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('chrom'),
        value=chromosomes,
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('sampleinfo', 'cell_id', axes_origin=[]),
        value=sample_info)


    workflow.setobj(
        obj=mgd.OutputChunks('cell_id', 'lane'),
        value=fastq_1_filename.keys(),
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('trim', 'cell_id', 'lane', axes_origin=[]),
        value=triminfo)

    workflow.setobj(
        obj=mgd.TempOutputObj('center', 'cell_id', 'lane', axes_origin=[]),
        value=centerinfo)

    fastqc_reports = os.path.join(
        lane_metrics,
        "fastqc",
        "{cell_id}_reports.tar.gz")
    flagstat_metrics = os.path.join(lane_metrics, 'flagstat', '{cell_id}.txt')
    workflow.transform(
        name='align_reads',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        axes=('cell_id', 'lane',),
        func="single_cell.workflows.align.tasks.align_pe",
        args=(
            mgd.InputFile(
                'fastq_1', 'cell_id', 'lane', fnames=fastq_1_filename),
            mgd.InputFile(
                'fastq_2', 'cell_id', 'lane', fnames=fastq_2_filename),
            mgd.TempOutputFile('biobloom_count_metrics', 'cell_id', 'lane'),
            mgd.TempSpace('alignment_temp', 'cell_id', 'lane'),
            ref_genome,
            config['biobloom_filters'],
        )
    )

    workflow.transform(
        name='merge_biobloom',
        func="single_cell.workflows.align.tasks.merge_biobloom",
        axes=('cell_id',),
        args=( mgd.TempInputFile('biobloom_count_metrics', 'cell_id', 'lane'),
               mgd.TempOutputFile('biobloom_count_metrics_merged', 'cell_id'),
               disable_biobloom,
               )
    )


    return workflow
