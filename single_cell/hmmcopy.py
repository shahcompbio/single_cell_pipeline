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

    for params_tag, params in config["hmmcopy_params"].iteritems():
        params_tag = "hmmcopy_" + params_tag

        results_dir = os.path.join(args['out_dir'], 'results', params_tag)
        plots_dir = os.path.join(results_dir, "plots")

        info_file = os.path.join(results_dir, "info.yaml")


        igv_seg_file = os.path.join(results_dir, '{}_igv_segments.seg'.format(lib))

        hmmcopy_data = os.path.join(results_dir, '{}_hmmcopy.h5'.format(lib))

        segs_pdf = os.path.join(plots_dir, lib+'_segs_row_{row}.pdf')
        bias_pdf = os.path.join(plots_dir, lib+'_bias_row_{row}.pdf')

        heatmap_filt_pdf = os.path.join(plots_dir, '{}_heatmap_by_ec_filtered.pdf'.format(lib))
        heatmap_pdf = os.path.join(plots_dir, '{}_heatmap_by_ec.pdf'.format(lib))
        metrics_pdf = os.path.join(plots_dir, '{}_metrics.pdf'.format(lib))
        kernel_density_pdf = os.path.join(plots_dir, '{}_kernel_density.pdf'.format(lib))

        workflow.subworkflow(
            name='hmmcopy_workflow_' + params_tag,
            func=hmmcopy.create_hmmcopy_workflow,
            args=(
                mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_files),
                mgd.InputFile('bam_markdups_index', 'cell_id', fnames=bai_files),
                mgd.OutputFile(hmmcopy_data),
                mgd.OutputFile(igv_seg_file),
                segs_pdf,
                bias_pdf,
                mgd.OutputFile(heatmap_pdf),
                mgd.OutputFile(heatmap_filt_pdf),
                mgd.OutputFile(metrics_pdf),
                mgd.OutputFile(kernel_density_pdf),
                mgd.OutputFile(info_file),
                cellids,
                config,
                args,
                params,
                params_tag,
                results_dir
            ),
        )

    return workflow
