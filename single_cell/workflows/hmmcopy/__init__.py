'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd

import single_cell.utils
import tasks


def create_hmmcopy_workflow(bam_file, bai_file, alignment_metrics, corrected_reads_file,
                            segments_file, reads_pdf_output, segs_pdf_output, bias_pdf_output,
                            params_pdf_output, hmm_metrics_file, sample_info,
                            sample_ids, config, args, hmmparams, results_dir ):

    lib = args['library_id']
    reads_filt_filename = os.path.join(results_dir, '{}_filtered_reads.csv'.format(lib))
    segs_filt_filename = os.path.join(results_dir, '{}_filtered_segs.csv'.format(lib))
    cn_matrix_file = os.path.join(results_dir, '{}_cn_matrix.csv'.format(lib))
    output_seg_filename = os.path.join(results_dir, '{}_igv_segments.seg'.format(lib))

    hmm_params_file = os.path.join(results_dir, '{}_params.csv'.format(lib))
    hmm_posteriors_file = os.path.join(results_dir, '{}_posteriors.csv'.format(lib))

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )

    workflow.transform(
        name='run_hmmcopy',
        ctx={'mem': config["memory"]['low'], 'pool_id': config['pools']['standard']},
        func=tasks.run_hmmcopy,
        axes=('sample_id',),
        args=(
            mgd.InputFile('bam_markdups', 'sample_id', fnames=bam_file),
            mgd.InputFile('bai_markdups', 'sample_id', fnames=bai_file),
            mgd.TempOutputFile('reads.csv', 'sample_id'),
            mgd.TempOutputFile('segs.csv', 'sample_id'),
            mgd.TempOutputFile('params.csv', 'sample_id'),
            mgd.TempOutputFile('posteriors.csv', 'sample_id'),
            mgd.TempOutputFile('hmm_metrics.csv', 'sample_id'),
            mgd.InputInstance('sample_id'),
            config,
            hmmparams,
            mgd.TempSpace('hmmcopy_temp', 'sample_id')
        ),
    )


    workflow.transform(
        name='merge_files',
        ctx={'mem': config["memory"]['low'], 'pool_id': config['pools']['standard']},
        func=tasks.merge_files,
        args=(
            mgd.TempInputFile('reads.csv', 'sample_id'),
            mgd.TempInputFile('segs.csv', 'sample_id'),
            mgd.TempInputFile('hmm_metrics.csv', 'sample_id'),
            mgd.TempInputFile('params.csv', 'sample_id'),
            mgd.TempInputFile('posteriors.csv', 'sample_id'),
            mgd.OutputFile(segments_file),
            mgd.OutputFile(corrected_reads_file),
            mgd.OutputFile(hmm_metrics_file),
            mgd.OutputFile(hmm_params_file),
            mgd.OutputFile(hmm_posteriors_file),
            mgd.OutputFile(cn_matrix_file),
            mgd.TempSpace('cnmatrix_temp'),
            0.2,
            mgd.OutputFile(reads_filt_filename),
            mgd.OutputFile(segs_filt_filename),
            mgd.OutputFile(output_seg_filename),
        ),
    )

    workflow.transform(
        name = 'merge_metrics',
        ctx={'mem': config["memory"]['low'], 'pool_id': config['pools']['standard']},
        func=tasks.merge_csv,
        args=(
             [mgd.InputFile(hmm_metrics_file),
              mgd.InputFile(alignment_metrics),
              mgd.InputFile(sample_info)],
             mgd.TempOutputFile("merged_metrics.csv"),
             'outer',
             'cell_id'
        )
    )

    workflow.transform(
        name='plot_hmm_copy',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard']},
        func=tasks.plot_hmmcopy,
        axes=('sample_id',),
        args=(
            mgd.TempInputFile('reads.csv', 'sample_id'),
            mgd.TempInputFile('segs.csv', 'sample_id'),
            mgd.TempInputFile('params.csv', 'sample_id'),
            mgd.TempInputFile("merged_metrics.csv"),
            mgd.InputFile(sample_info),
            config['ref_genome'],
            mgd.TempOutputFile('reads.pdf', 'sample_id'),
            mgd.TempOutputFile('segs.pdf', 'sample_id'),
            mgd.TempOutputFile('bias.pdf', 'sample_id'),
            mgd.TempOutputFile('params.pdf', 'sample_id'),
            mgd.InputInstance('sample_id'),
        ),
        kwargs={
            'num_states': hmmparams['num_states'],
            'plot_title': 'QC pipeline metrics',
            'annotation_cols': ['cell_call','experimental_condition', 'sample_type',
                                'coverage_depth', 'mad_neutral_state', 'MSRSI_non_integerness'],
        }
    )

    workflow.transform(
        name='merge_hmm_copy',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard']},
        func=tasks.merge_pdf,
        args=(
              [mgd.TempInputFile('reads.pdf', 'sample_id'),
              mgd.TempInputFile('segs.pdf', 'sample_id'),
              mgd.TempInputFile('bias.pdf', 'sample_id'),
              mgd.TempInputFile('params.pdf', 'sample_id')],
              [mgd.OutputFile(reads_pdf_output),
              mgd.OutputFile(segs_pdf_output),
              mgd.OutputFile(bias_pdf_output),
              mgd.OutputFile(params_pdf_output)],
              mgd.TempInputFile("merged_metrics.csv"),
              None
            )
    )

    reads_mad_pdf_output = os.path.join(results_dir, 'plots', '{}_reads_mad.pdf'.format(lib))
    segs_mad_pdf_output = os.path.join(results_dir, 'plots', '{}_segs_mad.pdf'.format(lib))
    bias_mad_pdf_output = os.path.join(results_dir, 'plots', '{}_bias_mad.pdf'.format(lib))
    params_mad_pdf_output = os.path.join(results_dir, 'plots', '{}_params_mad.pdf'.format(lib))
    workflow.transform(
        name='merge_hmm_copy_mad',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard']},
        func=tasks.merge_pdf,
        args=(
              [mgd.TempInputFile('reads.pdf', 'sample_id'),
              mgd.TempInputFile('segs.pdf', 'sample_id'),
              mgd.TempInputFile('bias.pdf', 'sample_id'),
              mgd.TempInputFile('params.pdf', 'sample_id')],
              [mgd.OutputFile(reads_mad_pdf_output),
              mgd.OutputFile(segs_mad_pdf_output),
              mgd.OutputFile(bias_mad_pdf_output),
              mgd.OutputFile(params_mad_pdf_output)],
              mgd.TempInputFile("merged_metrics.csv"),
              0.2
            )
    )


    return workflow
