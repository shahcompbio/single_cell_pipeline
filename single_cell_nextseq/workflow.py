

def create_alignment_workflow(
        fastq_1_filename,
        fastq_2_filename,
        bam_filename,
        bam_index_filename,
        ref_genome,
        read_group,
        config):

    workflow = pypeliner.workflow.Workflow()

    workflow.commandline(
        name='aln_read_1',
        ctx={'mem': 6},
        args=(
            'bwa',
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
            'bwa', 'aln',
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
            'bwa', 'sampe',
            '-r', read_group,
            mgd.InputFile(ref_genome),
            mgd.TempInputFile('read_1.sai'),
            mgd.TempInputFile('read_2.sai'),
            mgd.InputFile(fastq_1_filename),
            mgd.InputFile(fastq_2_filename),
            '|',
            'samtools', 'view',
            '-bSh', '-',
            '>',
            mgd.TempOutputFile('aligned.bam'),
        ),
    )

    workflow.transform(
        name='bam_sort',
        ctx={'mem': 4},
        func=tasks.bam_sort,
        args=(
            mgd.TempInputFile('aligned.bam'),
            mgd.TempOutputFile('sorted.bam'),
        ),
    )

    workflow.transform(
        name='bam_mark_dups',
        ctx={'mem': 4},
        func=tasks.bam_mark_dups,
        args=(
            mgd.TempInputFile('sorted.bam'),
            mgd.OutputFile(bam_filename),
            mgd.TempOutputFile('metrics.txt'),
        ),
    )

    workflow.commandline(
        name='bam_index',
        ctx={'mem': 4},
        args=(
            'samtools', 'index',
            mgd.InputFile(bam_filename),
            mgd.OutputFile(bam_index_filename),
        ),
    )

    return workflow


def create_metrics_workflow(
    bam_filename,
    ref_genome,
    metrics_summary_filename,
    metrics_directory):

    flagstat_metrics_filename = os.path.join(metrics_directory, 'flagstat_metrics/\2.flagstat_metrics.txt')
    wgs_metrics_filename = os.path.join(metrics_directory, 'wgs_metrics/\2.wgs_metrics.txt')
    gc_metrics_filename = os.path.join(metrics_directory, 'gc_metrics/\2.gc_metrics.txt')
    gc_summary_filename = os.path.join(metrics_directory, 'gc_metrics/\2.gc_metrics.summ.txt')
    gc_chart_filename = os.path.join(metrics_directory, 'gc_metrics/\2.gc_metrics.pdf')
    insert_metrics_filename = os.path.join(metrics_directory, 'insert_metrics/\2.insert_metrics.txt')
    insert_histogram_filename = os.path.join(metrics_directory, 'insert_metrics/\2.insert_metrics.pdf')

    workflow = pypeliner.workflow.Workflow()

    workflow.commandline(
        name='bam_flagstat',
        ctx={'mem': 4},
        args=(
            'samtools', 'flagstat',
            mgd.InputFile(bam_filename),
            '>',
            mgd.OutputFile(flagstat_metrics_filename),
        ),
    )


    workflow.commandline(
        name='bam_collect_wgs_metrics',
        ctx={'mem': 4},
        func=tasks.bam_collect_wgs_metrics,
        args=(
            mgd.InputFile(bam_filename),
            mgd.InputFile(ref_genome),
            mgd.OutputFile(wgs_metrics_filename),
            config,
        ),
    )

    workflow.commandline(
        name='bam_collect_gc_metrics',
        ctx={'mem': 4},
        func=tasks.bam_collect_gc_metrics,
        args=(
            mgd.InputFile(bam_filename),
            mgd.InputFile(ref_genome),
            mgd.OutputFile(gc_metrics_filename),
            mgd.OutputFile(gc_summary_filename),
            mgd.OutputFile(gc_chart_filename),
        ),
    )

    workflow.commandline(
        name='bam_collect_insert_metrics',
        ctx={'mem': 4},
        func=tasks.bam_collect_insert_metrics,
        args=(
            mgd.InputFile(bam_filename),
            mgd.OutputFile(insert_metrics_filename),
            mgd.OutputFile(insert_histogram_filename),
        ),
    )

    return workflow


