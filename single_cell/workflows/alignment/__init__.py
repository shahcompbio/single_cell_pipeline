'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import tasks
import pypeliner
import pypeliner.managed as mgd


def create_alignment_workflow(
        fastq_1_filename,
        fastq_2_filename,
        bam_filename,
        bai_filename,
        alignment_metrics,
        gc_metrics,
        ref_genome,
        config,
        args,
        seqinfo,
        sample_ids):

    out_dir = args['out_dir']

    merge_metrics = os.path.join(out_dir, 'metrics')
    lane_metrics = os.path.join(args['out_dir'], 'metrics_per_lane', '{lane}')

    bam_filename = dict([(sampid, bam_filename[sampid])
                         for sampid in sample_ids])

    bai_filename = dict([(sampid, bai_filename[sampid])
                         for sampid in sample_ids])

    chromosomes = config["chromosomes"]

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('chrom'),
        value=chromosomes,
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'lane'),
        value=fastq_1_filename.keys(),
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('seqinfo', 'sample_id', axes_origin=[]),
        value=seqinfo)


    fastqc_reports = os.path.join(lane_metrics, "fastqc", "{sample_id}_reports.tar.gz")
    flagstat_metrics = os.path.join(lane_metrics, 'flagstat', '{sample_id}.txt')
    workflow.transform(
        name='align_reads',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard']},
        axes = ('sample_id', 'lane',),
        func=tasks.align_pe,
        args=(
            mgd.InputFile('fastq_1', 'sample_id', 'lane', fnames=fastq_1_filename),
            mgd.InputFile('fastq_2', 'sample_id', 'lane', fnames=fastq_2_filename),
            mgd.TempOutputFile('aligned_per_cell_per_lane.sorted.bam', 'sample_id', 'lane'),
            mgd.OutputFile(fastqc_reports, 'sample_id', 'lane'),
            mgd.OutputFile(flagstat_metrics, 'sample_id', 'lane'),
            mgd.TempSpace('alignment_temp', 'sample_id', 'lane'),
            ref_genome,
            mgd.TempInputObj('seqinfo', 'sample_id'),
            mgd.InputInstance('sample_id'),
            mgd.InputInstance('lane'),
            args['library_id'],
            config
        )
    )

    workflow.transform(
        name='merge_bams',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard']},
        func=tasks.merge_bams,
        axes=('sample_id',),
        args=(
            mgd.TempInputFile('aligned_per_cell_per_lane.sorted.bam', 'sample_id', 'lane'),
            mgd.TempOutputFile('merged_lanes.bam', 'sample_id'),
            mgd.TempOutputFile('merged_lanes.bam.bai', 'sample_id'),
        )
    )

    if args['realign']:
        workflow.transform(
            name='realignment',
            axes=('chrom',),
            ctx={'mem': config["memory"]['high'], 'pool_id': config['pools']['highmem']},
            func=tasks.realign,
            args=(
                mgd.TempInputFile('merged_lanes.bam', 'sample_id'),
                mgd.TempInputFile('merged_lanes.bam.bai', 'sample_id'),
                mgd.TempOutputFile('realigned.bam', 'chrom', 'sample_id'),
                mgd.TempSpace('realignment_temp', 'chrom', cleanup='before'),
                config,
                mgd.InputInstance('chrom')
            )
        )

        workflow.transform(
            name='merge_realignment',
            ctx={'mem': config["memory"]['high'], 'pool_id': config['pools']['highmem']},
            axes=('sample_id',),
            func=tasks.merge_realignment,
            args=(
                mgd.TempInputFile('realigned.bam', 'chrom', 'sample_id'),
                mgd.TempOutputFile('merged_realign.bam', 'sample_id'),
                config,
                mgd.InputInstance('sample_id')
            )
        )

    final_bam = mgd.TempInputFile('merged_lanes.bam', 'sample_id')
    if args["realign"]:
        final_bam = mgd.TempInputFile('merged_realign.bam', 'sample_id')

    markdups_metrics = os.path.join(merge_metrics, 'markdups_metrics', '{sample_id}.markdups_metrics.txt')
    flagstat_metrics = os.path.join(merge_metrics, 'flagstat_metrics', '{sample_id}.flagstat_metrics.txt')
    workflow.transform(
        name='postprocess_bam',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard']},
        axes=('sample_id',),
        func=tasks.postprocess_bam,
        args=(
              final_bam,
            mgd.OutputFile('sorted_markdups', 'sample_id', fnames=bam_filename),
            mgd.OutputFile('sorted_markdups_index', 'sample_id', fnames=bai_filename),
            mgd.TempSpace('tempdir', 'sample_id'),
            config,
            mgd.OutputFile(markdups_metrics, 'sample_id'),
            mgd.OutputFile(flagstat_metrics, 'sample_id'),
        ),
    )

    wgs_metrics_filename = os.path.join(merge_metrics, 'wgs_metrics', '{sample_id}.wgs_metrics.txt') 
    workflow.transform(
        name='bam_collect_wgs_metrics',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard']},
        func=tasks.bam_collect_wgs_metrics,
        axes=('sample_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'sample_id', fnames=bam_filename),
            ref_genome,
            mgd.OutputFile(wgs_metrics_filename, 'sample_id'),
            config,
            mgd.TempSpace('wgs_tempdir', 'sample_id'),
        ),
    )

    gc_metrics_filename = os.path.join(merge_metrics, 'gc_metrics', '{sample_id}.gc_metrics.txt')
    gc_summary_filename = os.path.join(merge_metrics, 'gc_metrics', '{sample_id}.gc_metrics.summ.txt')
    gc_chart_filename = os.path.join(merge_metrics, 'gc_metrics', '{sample_id}.gc_metrics.pdf')
    workflow.transform(
        name='bam_collect_gc_metrics',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard']},
        func=tasks.bam_collect_gc_metrics,
        axes=('sample_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'sample_id', fnames=bam_filename),
            ref_genome,
            mgd.OutputFile(gc_metrics_filename, 'sample_id'),
            mgd.OutputFile(gc_summary_filename, 'sample_id'),
            mgd.OutputFile(gc_chart_filename, 'sample_id'),
            mgd.TempSpace('gc_tempdir', 'sample_id'),
        ),
    )

    insert_metrics_filename = os.path.join(merge_metrics, 'insert_metrics', '{sample_id}.insert_metrics.txt')
    insert_histogram_filename = os.path.join(merge_metrics, 'insert_metrics', '{sample_id}.insert_metrics.pdf')
    workflow.transform(
        name='bam_collect_insert_metrics',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard']},
        func=tasks.bam_collect_insert_metrics,
        axes=('sample_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'sample_id', fnames=bam_filename),
            mgd.InputFile(flagstat_metrics, 'sample_id'),
            mgd.OutputFile(insert_metrics_filename, 'sample_id'),
            mgd.OutputFile(insert_histogram_filename, 'sample_id'),
            mgd.TempSpace('insert_tempdir', 'sample_id'),
        ),
    )

    workflow.transform(
        name="collect_gc_metrics",
        func=tasks.collect_gc,
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard']},
        args = (
            mgd.InputFile(gc_metrics_filename, 'sample_id', axes_origin=[]),
            mgd.OutputFile(gc_metrics),
            mgd.TempSpace("temp_gc")
        )
    )

    workflow.transform(
        name='collect_metrics',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard']},
        func=tasks.collect_metrics,
        args=(
            mgd.InputFile(flagstat_metrics, 'sample_id', axes_origin=[]),
            mgd.InputFile(markdups_metrics, 'sample_id', axes_origin=[]),
            mgd.InputFile(insert_metrics_filename, 'sample_id', axes_origin=[]),
            mgd.InputFile(wgs_metrics_filename, 'sample_id', axes_origin=[]),
            mgd.TempSpace("tempdir_collect_metrics"),
            mgd.OutputFile(alignment_metrics),
        ),
    )

    return workflow
