'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd

import tasks


def create_hmmcopy_workflow(bam_file, bai_file, reads_file,
                            segments_file, metrics_file, params_file,
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
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
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
        ctx={'mem': config["memory"]['low'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.merge_files,
        args=(
            mgd.TempInputFile('reads.csv', 'sample_id'),
            mgd.TempInputFile('segs.csv', 'sample_id'),
            mgd.TempInputFile('hmm_metrics.csv', 'sample_id'),
            mgd.TempInputFile('params.csv', 'sample_id'),
            mgd.TempInputFile('posteriors.csv', 'sample_id'),
            mgd.OutputFile(segments_file),
            mgd.OutputFile(reads_file),
            mgd.OutputFile(metrics_file),
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

    return workflow
