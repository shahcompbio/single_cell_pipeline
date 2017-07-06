'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd

import tasks


def create_hmmcopy_workflow(bam_file, corrected_reads_file,
                            segments_file, hmm_metrics_file,
                            cnmatrix_file, sample_id, config, args):

    scripts_dir = os.path.join(os.path.realpath(os.path.dirname(__file__)),
                               'scripts')

    extract_quality_metrics_script = os.path.join(scripts_dir,
                                                  'extract_quality_metrics.py')
    cn_metrics_script = os.path.join(scripts_dir, 'gen_cn_matrix.py')

    chromosomes = config['chromosomes']

    hmmcopy_directory = os.path.join(args['out_dir'], 'hmmcopy')
    posterior_marginals_filename = os.path.join(hmmcopy_directory,
                                                sample_id + '_posteriors.csv')

    params_filename = os.path.join(hmmcopy_directory,
                                   sample_id + '_parameters.csv')

    wig_file =  os.path.join(hmmcopy_directory,
                             sample_id + '_readcounter.wig')



    workflow = pypeliner.workflow.Workflow()

    workflow.commandline(
        name='count_reads',
        ctx={'mem': 4},
        args=(
            config['readcounter'],
            '-w', str(config['bin_size']),
            '-q', str(config['min_mqual']),
            '-c', ','.join(chromosomes),
            mgd.InputFile(bam_file),
            '>',
            mgd.OutputFile(wig_file),
        ),
    )

    workflow.transform(
        name='run_hmmcopy',
        ctx={'mem': 4},
        func=tasks.run_hmmcopy,
        args=(
            mgd.InputFile(wig_file),
            mgd.OutputFile(corrected_reads_file),
            mgd.OutputFile(segments_file),
            mgd.OutputFile(params_filename),
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
            '--hmmcopy_params', mgd.InputFile(params_filename),
            '--hmmcopy_corrected_reads', mgd.InputFile(
                corrected_reads_file),
            '--hmmcopy_segments', mgd.InputFile(segments_file),
            '--out_file', mgd.OutputFile(hmm_metrics_file),
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
            '--input', mgd.InputFile(corrected_reads_file),
            '--output', mgd.OutputFile(cnmatrix_file),
            '--sample_id', sample_id,
            '--type', 'hmmcopy_corrected_reads',
            '--column_name', 'integer_copy_number'
        ),
    )

    return workflow
