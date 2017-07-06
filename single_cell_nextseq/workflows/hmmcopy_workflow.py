'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd

import single_cell_nextseq.tasks




def create_hmmcopy_workflow(
    bam_filename,
    wig_filename,
    corrected_reads_filename,
    segments_filename,
    parameters_filename,
    hmm_metrics_filename,
    cnmatrix_filename,
    sample_id,
    config,
    args):

    scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
    extract_quality_metrics_script = os.path.join(scripts_directory, 'extract_quality_metrics.py')
    cn_metrics_script = os.path.join(scripts_directory, 'gen_cn_matrix.py')

    chromosomes = config['chromosomes']

    hmmcopy_directory = os.path.join(args['out_dir'], 'hmmcopy')
    posterior_marginals_filename = os.path.join(hmmcopy_directory, sample_id+'_posteriors.csv')


    workflow = pypeliner.workflow.Workflow()

    workflow.commandline(
        name='count_reads',
        ctx={'mem': 4},
        args=(
            config['readcounter'],
            '-w', str(config['bin_size']),
            '-q', str(config['min_mqual']),
            '-c', ','.join(chromosomes),
            mgd.InputFile(bam_filename),
            '>',
            mgd.OutputFile(wig_filename),
        ),
    )

    workflow.transform(
        name='run_hmmcopy',
        ctx={'mem': 4},
        func=single_cell_nextseq.tasks.run_hmmcopy,
        args=(
            mgd.InputFile(wig_filename),
            mgd.OutputFile(corrected_reads_filename),
            mgd.OutputFile(segments_filename),
            mgd.OutputFile(parameters_filename),
            mgd.OutputFile(posterior_marginals_filename),
            sample_id,
            config,
        ),
    )

    workflow.commandline(
        name='extract_quality_metrics',
        ctx={'mem': 4},
        args=(
            config['python'],
            extract_quality_metrics_script,
            '--hmmcopy_params', mgd.InputFile(parameters_filename),
            '--hmmcopy_corrected_reads', mgd.InputFile(corrected_reads_filename),
            '--hmmcopy_segments', mgd.InputFile(segments_filename),
            '--out_file', mgd.OutputFile(hmm_metrics_filename),
            '--sample_id', sample_id,
        ),
    )


    workflow.commandline(
        name='collect_gc_metrics',
        ctx={'mem': 16},
        args=(
            config['python'],
            cn_metrics_script,
            '--separator', 'comma',
            '--input', mgd.InputFile(corrected_reads_filename),
            '--output', mgd.OutputFile(cnmatrix_filename),
            '--sample_id', sample_id,
            '--type', 'hmmcopy_corrected_reads', 
            '--column_name', 'integer_copy_number'
        ),
    )

    return workflow

