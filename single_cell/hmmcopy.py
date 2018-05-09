'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import pypeliner.managed as mgd
from workflows import hmmcopy
from utils import helpers


def hmmcopy_workflow(workflow, args):

    config = helpers.load_config(args)
    cellids = helpers.get_samples(args['input_yaml'])
    bam_files, bai_files = helpers.get_bams(args['input_yaml'])
    lib = args['library_id']

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cellids,
    )

    for name, params in config["hmmcopy_params"].iteritems():
        name = "hmmcopy_" + name

        results_dir = os.path.join(args['out_dir'], 'results', name)
        plots_dir = os.path.join(results_dir, "plots")

        hmmcopy_data = os.path.join(results_dir, '{}_hmmcopy.h5'.format(lib))
        igv_seg_file = os.path.join(results_dir, '{}_igv_segments.seg'.format(lib))

        reads_pdf = os.path.join(plots_dir, '{}_reads.pdf'.format(lib))
        segs_pdf = os.path.join(plots_dir, '{}_segs.pdf'.format(lib))
        bias_pdf = os.path.join(plots_dir, '{}_bias.pdf'.format(lib))
        params_pdf = os.path.join(plots_dir, '{}_params.pdf'.format(lib))

        reads_filt_pdf = os.path.join(plots_dir, '{}_reads_filtered.pdf'.format(lib))
        segs_filt_pdf = os.path.join(plots_dir, '{}_segs_filtered.pdf'.format(lib))
        bias_filt_pdf = os.path.join(plots_dir, '{}_bias_filtered.pdf'.format(lib))
        params_filt_pdf = os.path.join(plots_dir, '{}_params_filtered.pdf'.format(lib))

        heatmap_filt_pdf = os.path.join(plots_dir, '{}_heatmap_by_ec_filtered.pdf'.format(lib))
        heatmap_pdf = os.path.join(plots_dir, '{}_heatmap_by_ec.pdf'.format(lib))
        metrics_pdf = os.path.join(plots_dir, '{}_metrics.pdf'.format(lib))
        kernel_density_pdf = os.path.join(plots_dir, '{}_kernel_density.pdf'.format(lib))

        workflow.subworkflow(
            name='hmmcopy_workflow_' + name,
            func=hmmcopy.create_hmmcopy_workflow,
            args=(
                mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_files),
                mgd.InputFile('bam_markdups_index', 'cell_id', fnames=bai_files),
                mgd.OutputFile(hmmcopy_data),
                mgd.OutputFile(igv_seg_file),
                mgd.OutputFile(reads_pdf),
                mgd.OutputFile(segs_pdf),
                mgd.OutputFile(bias_pdf),
                mgd.OutputFile(params_pdf),
                mgd.OutputFile(reads_filt_pdf),
                mgd.OutputFile(segs_filt_pdf),
                mgd.OutputFile(bias_filt_pdf),
                mgd.OutputFile(params_filt_pdf),
                mgd.OutputFile(heatmap_pdf),
                mgd.OutputFile(heatmap_filt_pdf),
                mgd.OutputFile(metrics_pdf),
                mgd.OutputFile(kernel_density_pdf),
                cellids,
                config,
                args,
                params,
                results_dir
            ),
        )

    return workflow
