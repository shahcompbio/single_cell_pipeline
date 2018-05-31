'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd


def create_alignment_workflow(
        fastq_1_filename,
        fastq_2_filename,
        bam_filename,
        bai_filename,
        alignment_metrics,
        plot_metrics,
        meta_yaml,
        ref_genome,
        config,
        args,
        instrumentinfo,
        centerinfo,
        sample_info,
        cell_ids):

    out_dir = args['out_dir']

    merge_metrics = os.path.join(out_dir, 'metrics')

    lane_metrics = os.path.join(args['out_dir'], 'metrics_per_lane', '{lane}')

    bam_filename = dict([(cellid, bam_filename[cellid])
                         for cellid in cell_ids])

    bai_filename = dict([(cellid, bai_filename[cellid])
                         for cellid in cell_ids])

    chromosomes = config["chromosomes"]

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('chrom'),
        value=chromosomes,
    )

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id', 'lane'),
        value=fastq_1_filename.keys(),
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('instrument', 'cell_id', 'lane', axes_origin=[]),
        value=instrumentinfo)

    workflow.setobj(
        obj=mgd.TempOutputObj('center', 'cell_id', 'lane', axes_origin=[]),
        value=centerinfo)


    workflow.setobj(
        obj=mgd.TempOutputObj('sampleinfo', 'cell_id', axes_origin=[]),
        value=sample_info)

    fastqc_reports = os.path.join(
        lane_metrics,
        "fastqc",
        "{cell_id}_reports.tar.gz")
    flagstat_metrics = os.path.join(lane_metrics, 'flagstat', '{cell_id}.txt')
    workflow.transform(
        name='align_reads',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        axes=('cell_id', 'lane',),
        func="single_cell.workflows.alignment.tasks.align_pe",
        args=(
            mgd.InputFile(
                'fastq_1',
                'cell_id',
                'lane',
                fnames=fastq_1_filename),
            mgd.InputFile(
                'fastq_2',
                'cell_id',
                'lane',
                fnames=fastq_2_filename),
            mgd.TempOutputFile(
                'aligned_per_cell_per_lane.sorted.bam',
                'cell_id',
                'lane'),
            mgd.OutputFile(fastqc_reports, 'cell_id', 'lane'),
            mgd.OutputFile(flagstat_metrics, 'cell_id', 'lane'),
            mgd.TempSpace('alignment_temp', 'cell_id', 'lane'),
            ref_genome,
            mgd.TempInputObj('instrument', 'cell_id', 'lane'),
            mgd.TempInputObj('center', 'cell_id', 'lane'),
            mgd.TempInputObj('sampleinfo', 'cell_id'),
            mgd.InputInstance('cell_id'),
            mgd.InputInstance('lane'),
            args['library_id'],
            config
        )
    )

    workflow.transform(
        name='merge_bams',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.alignment.tasks.merge_bams",
        axes=('cell_id',),
        args=(
            mgd.TempInputFile(
                'aligned_per_cell_per_lane.sorted.bam',
                'cell_id',
                'lane'),
            mgd.TempOutputFile('merged_lanes.bam', 'cell_id'),
            mgd.TempOutputFile('merged_lanes.bam.bai', 'cell_id'),
        )
    )

    if args['realign']:
        workflow.transform(
            name='realignment',
            axes=('chrom',),
            ctx={
                'mem': config["memory"]['high'],
                'pool_id': config['pools']['highmem'],
                'ncpus': 1},
            func="single_cell.workflows.alignment.tasks.realign",
            args=(
                mgd.TempInputFile('merged_lanes.bam', 'cell_id'),
                mgd.TempInputFile('merged_lanes.bam.bai', 'cell_id'),
                mgd.TempOutputFile('realigned.bam', 'chrom', 'cell_id'),
                mgd.TempSpace('realignment_temp', 'chrom', cleanup='before'),
                config,
                mgd.InputInstance('chrom')
            )
        )

        workflow.transform(
            name='merge_realignment',
            ctx={
                'mem': config["memory"]['high'],
                'pool_id': config['pools']['highmem'],
                'ncpus': 1},
            axes=('cell_id',),
            func="single_cell.workflows.alignment.tasks.merge_realignment",
            args=(
                mgd.TempInputFile('realigned.bam', 'chrom', 'cell_id'),
                mgd.TempOutputFile('merged_realign.bam', 'cell_id'),
                config,
                mgd.InputInstance('cell_id')
            )
        )

    final_bam = mgd.TempInputFile('merged_lanes.bam', 'cell_id')
    if args["realign"]:
        final_bam = mgd.TempInputFile('merged_realign.bam', 'cell_id')

    markdups_metrics = os.path.join(
        merge_metrics,
        'markdups_metrics',
        '{cell_id}.markdups_metrics.txt')
    flagstat_metrics = os.path.join(
        merge_metrics,
        'flagstat_metrics',
        '{cell_id}.flagstat_metrics.txt')
    workflow.transform(
        name='postprocess_bam',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        axes=('cell_id',),
        func="single_cell.workflows.alignment.tasks.postprocess_bam",
        args=(
            final_bam,
            mgd.OutputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            mgd.OutputFile(
                'sorted_markdups_index',
                'cell_id',
                fnames=bai_filename),
            mgd.TempSpace('tempdir', 'cell_id'),
            config,
            mgd.OutputFile(markdups_metrics, 'cell_id'),
            mgd.OutputFile(flagstat_metrics, 'cell_id'),
        ),
    )

    wgs_metrics_filename = os.path.join(
        merge_metrics,
        'wgs_metrics',
        '{cell_id}.wgs_metrics.txt')
    workflow.transform(
        name='bam_collect_wgs_metrics',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.alignment.tasks.bam_collect_wgs_metrics",
        axes=('cell_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            ref_genome,
            mgd.OutputFile(wgs_metrics_filename, 'cell_id'),
            config,
            mgd.TempSpace('wgs_tempdir', 'cell_id'),
        ),
    )

    gc_metrics_filename = os.path.join(
        merge_metrics,
        'gc_metrics',
        '{cell_id}.gc_metrics.txt')
    gc_summary_filename = os.path.join(
        merge_metrics,
        'gc_metrics',
        '{cell_id}.gc_metrics.summ.txt')
    gc_chart_filename = os.path.join(
        merge_metrics,
        'gc_metrics',
        '{cell_id}.gc_metrics.pdf')
    workflow.transform(
        name='bam_collect_gc_metrics',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.alignment.tasks.bam_collect_gc_metrics",
        axes=('cell_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            ref_genome,
            mgd.OutputFile(gc_metrics_filename, 'cell_id'),
            mgd.OutputFile(gc_summary_filename, 'cell_id'),
            mgd.OutputFile(gc_chart_filename, 'cell_id'),
            mgd.TempSpace('gc_tempdir', 'cell_id'),
        ),
    )

    insert_metrics_filename = os.path.join(
        merge_metrics,
        'insert_metrics',
        '{cell_id}.insert_metrics.txt')
    insert_histogram_filename = os.path.join(
        merge_metrics,
        'insert_metrics',
        '{cell_id}.insert_metrics.pdf')
    workflow.transform(
        name='bam_collect_insert_metrics',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.alignment.tasks.bam_collect_insert_metrics",
        axes=('cell_id',),
        args=(
            mgd.InputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            mgd.InputFile(flagstat_metrics, 'cell_id'),
            mgd.OutputFile(insert_metrics_filename, 'cell_id'),
            mgd.OutputFile(insert_histogram_filename, 'cell_id'),
            mgd.TempSpace('insert_tempdir', 'cell_id'),
        ),
    )

    workflow.transform(
        name="collect_gc_metrics",
        func="single_cell.workflows.alignment.tasks.collect_gc",
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        args=(
            mgd.InputFile(gc_metrics_filename, 'cell_id', axes_origin=[]),
            mgd.TempOutputFile("gc_metrics.h5"),
            mgd.TempSpace("temp_gc")
        )
    )

    workflow.transform(
        name='collect_metrics',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.alignment.tasks.collect_metrics",
        args=(
            mgd.InputFile(flagstat_metrics, 'cell_id', axes_origin=[]),
            mgd.InputFile(markdups_metrics, 'cell_id', axes_origin=[]),
            mgd.InputFile(insert_metrics_filename, 'cell_id', axes_origin=[]),
            mgd.InputFile(wgs_metrics_filename, 'cell_id', axes_origin=[]),
            mgd.TempSpace("tempdir_collect_metrics"),
            mgd.TempOutputFile("alignment_metrics.h5"),
        ),
    )

    workflow.transform(
        name='annotate_metrics',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.alignment.tasks.annotate_metrics",
        args=(
            mgd.TempInputFile("alignment_metrics.h5"),
            sample_info,
            mgd.TempOutputFile("alignment_metrics_annotated.h5"),
        )
    )

    workflow.transform(
        name='plot_metrics',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.alignment.tasks.plot_metrics",
        args=(
            mgd.TempInputFile("alignment_metrics_annotated.h5"),
            mgd.OutputFile(plot_metrics),
            'QC pipeline metrics',
            mgd.TempInputFile("gc_metrics.h5"),
            config['gc_windows'],
        )
    )

    workflow.transform(
        name='concatenate_all_hdf_tables',
        ctx={
            'mem': config["memory"]['low'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.workflows.alignment.hdfutils.concat_hdf_tables",
        args=(
            [mgd.TempInputFile("alignment_metrics_annotated.h5"),
             mgd.TempInputFile("gc_metrics.h5"),
             ],
            mgd.OutputFile(alignment_metrics),
        ),
    )


    metadata = {
        'alignment':{
            'cell_batch_realign': args["realign"],
            'data': alignment_metrics,
            'metrics_table': '/alignment/metrics',
            'gc_metrics_table': '/alignment/gc_metrics',
            'plots': plot_metrics,
            'aligner': config["aligner"],
            'adapter': config["adapter"],
            'adapter2': config["adapter2"],
            'picardtools_wgsmetrics_params': config['picard_wgs_params'],
            'ref_genome': config["ref_genome"]
        }
    }


    workflow.transform(
        name='generate_meta_yaml',
        ctx={
            'mem': config["memory"]['low'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func="single_cell.utils.helpers.write_to_yaml",
        args=(
            mgd.OutputFile(meta_yaml),
            metadata
        )
    )

    return workflow
