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
                            cnmatrix_file, sample_id, config, out_dir):

    hmmcopy_directory = os.path.join(out_dir, 'hmmcopy', 'intermediates')

    posterior_marginals_filename = os.path.join(hmmcopy_directory, '{}_posteriors.csv'.format(sample_id))
 
    params_filename = os.path.join(hmmcopy_directory, '{}_parameters.csv'.format(sample_id))

    wig_file =  os.path.join(hmmcopy_directory, '{}_readcounter.wig'.format(sample_id))

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='run_hmmcopy',
        ctx={'mem': config['low_mem']},
        func=tasks.run_hmmcopy,
        args=(
            mgd.InputFile(wig_file),
            mgd.OutputFile(corrected_reads_file),
            mgd.OutputFile(segments_file),
            mgd.OutputFile(params_filename),
            mgd.OutputFile(posterior_marginals_filename),
            mgd.OutputFile(hmm_metrics_file),
            sample_id,
            config,
            mgd.TempSpace('hmmcopy_temp')
        ),
    )
 
    workflow.transform(
        name='collect_cn_metrics',
        ctx={'mem': config['low_mem']},
        func=tasks.generate_cn_matrix,
        args=(
            mgd.InputFile(corrected_reads_file),
            mgd.OutputFile(cnmatrix_file),
            ',',
            'integer_copy_number',
            sample_id,
            'hmmcopy_corrected_reads',
        ),
    )

    return workflow
