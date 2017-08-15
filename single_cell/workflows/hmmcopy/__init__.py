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
                            sample_ids, config, args):

    results_dir = os.path.join(args['out_dir'], 'results')
    reads_filt_filename = os.path.join(results_dir, 'filtered_reads.csv')
    segs_filt_filename = os.path.join(results_dir, 'filtered_segs.csv')
    cn_matrix_file = os.path.join(results_dir, 'cn_matrix.csv')


    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )

    workflow.transform(
        name='run_hmmcopy',
        ctx={'mem': config['low_mem']},
        func=tasks.run_hmmcopy,
        axes=('sample_id',),
        args=(
            mgd.InputFile('bam_markdups', 'sample_id', fnames=bam_file),
            mgd.TempOutputFile('reads.csv', 'sample_id'),
            mgd.TempOutputFile('segs.csv', 'sample_id'),
            mgd.TempOutputFile('params.csv', 'sample_id'),
            mgd.TempOutputFile('posteriors.csv', 'sample_id'),
            mgd.TempOutputFile('hmm_metrics.csv', 'sample_id'),
            mgd.InputInstance('sample_id'),
            config,
            mgd.TempSpace('hmmcopy_temp', 'sample_id')
        ),
    )


    workflow.transform(
        name='plot_hmm_copy',
        ctx={'mem': config['high_mem']},
        func=tasks.plot_hmmcopy,
        axes=('sample_id',),
        args=(
            mgd.TempInputFile('reads.csv', 'sample_id'),
            mgd.TempInputFile('segs.csv', 'sample_id'),
            mgd.TempInputFile('hmm_metrics.csv', 'sample_id'),
            mgd.InputFile(config['ref_genome']),
            mgd.TempOutputFile('reads.pdf', 'sample_id'),
            mgd.TempOutputFile('segs.pdf', 'sample_id'),
            mgd.TempOutputFile('bias.pdf', 'sample_id'),
            mgd.InputInstance('sample_id')
        ),
        kwargs={
            'num_states': config['hmmcopy_params']['num_states'],
            'plot_title': 'QC pipeline metrics',
        }
    )

    workflow.transform(
        name='generate_cn_matrix',
        ctx={'mem': config['low_mem']},
        func=tasks.generate_cn_matrix,
        args=(
            mgd.TempInputFile('reads.csv', 'sample_id'),
            mgd.OutputFile(cn_matrix_file),
            mgd.TempSpace('cnmatrix_temp')
        ),
    )

    workflow.transform(
        name='merge_tables',
        ctx={'mem': config['med_mem']},
        func=tasks.concatenate_csv,
        args=(
            mgd.TempInputFile('segs.csv', 'sample_id'),
            mgd.OutputFile(segments_file),
        ),
    )

    workflow.transform(
        name='merge_reads',
        ctx={'mem': config['high_mem']},
        func=tasks.concatenate_csv,
        args=(
            mgd.TempInputFile('reads.csv', 'sample_id'),
            mgd.OutputFile(corrected_reads_file),
        ),
    )

    workflow.transform(
        name='merge_hmm_metrics',
        ctx={'mem': config['low_mem']},
        func=tasks.concatenate_csv,
        args=(
            mgd.TempInputFile('hmm_metrics.csv', 'sample_id'),
            mgd.OutputFile(hmm_metrics_file),
        ),
    )

    workflow.transform(
        name='filter_hmmcopy_results',
        ctx={'mem': config['high_mem']},
        func=tasks.filter_hmm_data,
        args=(
            mgd.InputFile(hmm_metrics_file),
            mgd.InputFile(segments_file),
            mgd.InputFile(corrected_reads_file),
            0.2,
            mgd.OutputFile(reads_filt_filename),
            mgd.OutputFile(segs_filt_filename),
        )
    )

    reads_pdf_output = os.path.join(results_dir, 'plots', 'reads.pdf')
    segs_pdf_output = os.path.join(results_dir, 'plots', 'segs.pdf')
    bias_pdf_output = os.path.join(results_dir, 'plots', 'bias.pdf')
    workflow.transform(
        name='merge_hmm_copy',
        ctx={'mem': config['high_mem']},
        func=tasks.merge_pdf,
        args=(
              [mgd.TempInputFile('reads.pdf', 'sample_id'),
              mgd.TempInputFile('segs.pdf', 'sample_id'),
              mgd.TempInputFile('bias.pdf', 'sample_id')],
              [reads_pdf_output,
              segs_pdf_output,
              bias_pdf_output],
              mgd.InputFile(hmm_metrics_file),
              None
            )
    )


    reads_mad_pdf_output = os.path.join(results_dir, 'plots', 'reads_mad.pdf')
    segs_mad_pdf_output = os.path.join(results_dir, 'plots', 'segs_mad.pdf')
    bias_mad_pdf_output = os.path.join(results_dir, 'plots', 'bias_mad.pdf')
    workflow.transform(
        name='merge_hmm_copy_mad',
        ctx={'mem': config['high_mem']},
        func=tasks.merge_pdf,
        args=(
              [mgd.TempInputFile('reads.pdf', 'sample_id'),
              mgd.TempInputFile('segs.pdf', 'sample_id'),
              mgd.TempInputFile('bias.pdf', 'sample_id')],
              [mgd.OutputFile(reads_mad_pdf_output),
              mgd.OutputFile(segs_mad_pdf_output),
              mgd.OutputFile(bias_mad_pdf_output)],
              mgd.InputFile(hmm_metrics_file),
              0.2
            )
    )


    return workflow
