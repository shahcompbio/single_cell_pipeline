'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import pypeliner.managed as mgd
from workflows import copyclone
from utils import helpers


def copyclone_workflow(workflow, args):

    config = helpers.load_config(args)
    sampleids = helpers.get_samples(args['input_yaml'])
    bam_files, bai_files = helpers.get_bams(args['input_yaml'])

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sampleids,
    )

    libid = args['library_id']

    results_dir = os.path.join(args['out_dir'], 'results', "copyclone")
    output_filename = os.path.join(results_dir, '{}_copyclone.h5'.format(libid))

    plots_dir = os.path.join(results_dir, "plots")
    segs_pdf = os.path.join(plots_dir, '{}_segments.pdf'.format(libid))
    reads_pdf = os.path.join(plots_dir, '{}_bias.pdf'.format(libid))
    metrics_pdf = os.path.join(plots_dir, '{}_metrics.pdf'.format(libid))
    kde_pdf = os.path.join(plots_dir, '{}_kernel_density.pdf'.format(libid))
    heatmap_pdf = os.path.join(plots_dir, '{}_heatmap.pdf'.format(libid))
    heatmap_filt_pdf = os.path.join(
        plots_dir,
        '{}_heatmap_filtered.pdf'.format(libid))

    workflow.subworkflow(
        name='copyclone_workflow',
        func=copyclone.create_copyclone_workflow,
        args=(mgd.InputFile('bam_markdups', 'sample_id', fnames=bam_files),
              mgd.InputFile(
            'bam_markdups_index',
            'sample_id',
            fnames=bai_files),
            mgd.OutputFile(output_filename),
            mgd.OutputFile(segs_pdf),
            mgd.OutputFile(reads_pdf),
            mgd.OutputFile(metrics_pdf),
            mgd.OutputFile(kde_pdf),
            mgd.OutputFile(heatmap_pdf),
            mgd.OutputFile(heatmap_filt_pdf),
            sampleids,
            config,
            args,
            results_dir
        ),
    )

    return workflow
