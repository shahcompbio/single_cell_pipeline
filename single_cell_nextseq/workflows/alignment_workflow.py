'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd

import tasks

def create_alignment_workflow(
    fastq_1_filename,
    fastq_2_filename,
    bam_filename,
    bam_index_filename,
    ref_genome,
    read_group,
    metrics_summary_filename,
    gc_matrix_filename,
    metrics_directory,
    sample_id,
    samplesheet,
    config):
    

    scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
    collect_metrics_script = os.path.join(scripts_directory, 'collect_metrics.py')
    scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
    gc_metrics_script = os.path.join(scripts_directory, 'gen_cn_matrix.py')


    markdups_metrics_filename = os.path.join(metrics_directory, 'markdups_metrics/{}.markdups_metrics.txt'.format(sample_id))
    flagstat_metrics_filename = os.path.join(metrics_directory, 'flagstat_metrics/{}.flagstat_metrics.txt'.format(sample_id))
    wgs_metrics_filename = os.path.join(metrics_directory, 'wgs_metrics/{}.wgs_metrics.txt'.format(sample_id))
    gc_metrics_filename = os.path.join(metrics_directory, 'gc_metrics/{}.gc_metrics.txt'.format(sample_id))
    gc_summary_filename = os.path.join(metrics_directory, 'gc_metrics/{}.gc_metrics.summ.txt'.format(sample_id))
    gc_chart_filename = os.path.join(metrics_directory, 'gc_metrics/{}.gc_metrics.pdf'.format(sample_id))
    insert_metrics_filename = os.path.join(metrics_directory, 'insert_metrics/{}.insert_metrics.txt'.format(sample_id))
    insert_histogram_filename = os.path.join(metrics_directory, 'insert_metrics/{}.insert_metrics.pdf'.format(sample_id))

    workflow = pypeliner.workflow.Workflow()

    workflow.commandline(
        name='aln_read_1',
        ctx={'mem': 6},
        args=(
            config['bwa'],
            'aln',
            mgd.InputFile(ref_genome),
            mgd.InputFile(fastq_1_filename),
            '>',
            mgd.TempOutputFile('read_1.sai'),
        ),
    )

    workflow.commandline(
        name='aln_read_2',
        ctx={'mem': 6},
        args=(
            config['bwa'], 'aln',
            mgd.InputFile(ref_genome),
            mgd.InputFile(fastq_2_filename),
            '>',
            mgd.TempOutputFile('read_2.sai'),
        ),
    )

    workflow.commandline(
        name='sampe',
        ctx={'mem': 6},
        args=(
            config['bwa'], 'sampe',
            '-r', read_group,
            mgd.InputFile(ref_genome),
            mgd.TempInputFile('read_1.sai'),
            mgd.TempInputFile('read_2.sai'),
            mgd.InputFile(fastq_1_filename),
            mgd.InputFile(fastq_2_filename),
            '|',
            config['samtools'], 'view',
            '-bSh', '-',
            '>',
            mgd.TempOutputFile('aligned.bam'),
        ),
    )

    workflow.transform(
        name='bam_sort',
        ctx={'mem': 16},
        func=tasks.bam_sort,
        args=(
            mgd.TempInputFile('aligned.bam'),
            mgd.TempOutputFile('sorted.bam'),
            config
        ),
    )

    workflow.transform(
        name='bam_markdups',
        ctx={'mem': 16},
        func=tasks.bam_markdups,
        args=(
            mgd.TempInputFile('sorted.bam'),
            mgd.OutputFile(bam_filename),
            mgd.OutputFile(markdups_metrics_filename),
            config
        ),
    )

    workflow.commandline(
        name='bam_index',
        ctx={'mem': 16},
        args=(
            config['samtools'], 'index',
            mgd.InputFile(bam_filename),
            # mgd.OutputFile(bam_index_filename),
        ),
    )

    workflow.commandline(
        name='bam_flagstat',
        ctx={'mem': 16},
        args=(
            config['samtools'], 'flagstat',
            mgd.InputFile(bam_filename),
            '>',
            mgd.OutputFile(flagstat_metrics_filename),
        ),
    )

    workflow.transform(
        name='bam_collect_wgs_metrics',
        ctx={'mem': 24},
        func=tasks.bam_collect_wgs_metrics,
        args=(
            mgd.InputFile(bam_filename),
            mgd.InputFile(ref_genome),
            mgd.OutputFile(wgs_metrics_filename),
            config,
        ),
    )

    workflow.transform(
        name='bam_collect_gc_metrics',
        ctx={'mem': 24},
        func=tasks.bam_collect_gc_metrics,
        args=(
            mgd.InputFile(bam_filename),
            mgd.InputFile(ref_genome),
            mgd.OutputFile(gc_metrics_filename),
            mgd.OutputFile(gc_summary_filename),
            mgd.OutputFile(gc_chart_filename),
            config
        ),
    )

    workflow.transform(
        name='bam_collect_insert_metrics',
        ctx={'mem': 24},
        func=tasks.bam_collect_insert_metrics,
        args=(
            mgd.InputFile(bam_filename),
            mgd.InputFile(flagstat_metrics_filename),
            mgd.OutputFile(insert_metrics_filename),
            mgd.OutputFile(insert_histogram_filename),
            config
        ),
    )
    
    workflow.commandline(
        name='collect_metrics',
        ctx={'mem': 16},
        args=(
            config['python'],
            collect_metrics_script,
            mgd.InputFile(flagstat_metrics_filename),
            mgd.InputFile(markdups_metrics_filename),
            mgd.InputFile(insert_metrics_filename),
            mgd.InputFile(wgs_metrics_filename),
            mgd.OutputFile(metrics_summary_filename),
            mgd.InputFile(samplesheet),
            sample_id,
        ),
    )

    workflow.commandline(
        name='collect_gc_metrics',
        ctx={'mem': 16},
        args=(
            config['python'],
            gc_metrics_script,
            '--separator', 'comma',
            '--input', mgd.InputFile(gc_metrics_filename),
            '--output', mgd.OutputFile(gc_matrix_filename),
            '--sample_id', sample_id,
            '--type', 'gcbias', 
            '--column_name', 'NORMALIZED_COVERAGE'
        ),
    )


    return workflow
