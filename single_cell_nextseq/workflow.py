import pypeliner
import pypeliner.managed as mgd


scripts_directory = os.path.join(os.path.readlpath(os.dirname(__file__)), 'scripts')
collect_metrics_script = os.path.join(scripts_directory, 'collect_metrics.py')


def create_alignment_workflow(
    fastq_1_filename,
    fastq_2_filename,
    bam_filename,
    bam_index_filename,
    ref_genome,
    read_group,
    metrics_summary_filename,
    metrics_directory,
    config):

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
        name='bam_markdups',
        ctx={'mem': 4},
        func=tasks.bam_markdups,
        args=(
            mgd.TempInputFile('sorted.bam'),
            mgd.OutputFile(bam_filename),
            mgd.OutputFile(markdups_metrics_filename),
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
    
    workflow.commandline(
        name='collect_metrics',
        ctx={'mem': 4},
        args=(
            'python',
            collect_metrics_script,
            mgd.InputFile(flagstat_metrics_filename),
            mgd.InputFile(markdups_metrics_filename),
            mgd.InputFile(insert_metrics_filename),
            mgd.InputFile(wgs_metrics_filename),
            mgd.OutputFile(metrics_summary_filename),
        ),
    )

    return workflow


def create_hmmcopy_workflow(
    bam_filename,
    config):

    chromosomes = config['chromosomes']

    workflow = pypeliner.workflow.Workflow()

    workflow.commandline(
        name='count_reads',
        ctx={'mem': 4},
        args=(
            'readCounter',
            '-w', str(config['bin_size']),
            '-q', str(config['min_mqual']),
            '-c', chroms,
            mgd.InputFile(bam_filename),
            '>',
            mgd.OutputFile(wig_filename),
        ),
    )

    workflow.commandline(
        name='run_hmmcopy',
        ctx={'mem': 4},
        args=(
            mgd.InputFile(wig_filename),
            mgd.OutputFile(corrected_reads_filename),
            mgd.OutputFile(segments_filename),
            mgd.OutputFile(parameters_filename),
            mgd.OutputFile(posterior_marginals_filename),
            mgd.TempSpace('hmmcopy_temp'),
            sample_id,
            config,
        ),
    )

    return workflow

