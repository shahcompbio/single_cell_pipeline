'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import pypeliner.managed as mgd
from workflows import hmmcopy
from utils import helpers
import single_cell


def hmmcopy_workflow(workflow, args):

    config = helpers.load_config(args)

    ctx = {'mem_retry_increment': 2, 'ncpus': 1}
    ctx.update(helpers.get_container_ctx(config['containers'], 'single_cell_pipeline'))

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

        igv_seg_file = os.path.join(
            results_dir, '{}_igv_segments.seg'.format(lib))

        hmmcopy_data = os.path.join(results_dir, '{}_hmmcopy.h5'.format(lib))

        segs_pdf = os.path.join(
            plots_dir, "segments", lib + '_segs.tar.gz')
        bias_pdf = os.path.join(plots_dir, "bias", lib + '_bias.tar.gz')

        heatmap_filt_pdf = os.path.join(
            plots_dir, '{}_heatmap_by_ec_filtered.pdf'.format(lib))
        heatmap_pdf = os.path.join(
            plots_dir, '{}_heatmap_by_ec.pdf'.format(lib))
        metrics_pdf = os.path.join(plots_dir, '{}_metrics.pdf'.format(lib))
        kernel_density_pdf = os.path.join(
            plots_dir, '{}_kernel_density.pdf'.format(lib))

        workflow.subworkflow(
            name='hmmcopy_workflow_' + params_tag,
            func=hmmcopy.create_hmmcopy_workflow,
            args=(
                mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_files),
                mgd.InputFile(
                    'bam_markdups_index', 'cell_id', fnames=bai_files),
                mgd.OutputFile(hmmcopy_data),
                mgd.OutputFile(igv_seg_file),
                segs_pdf,
                bias_pdf,
                mgd.OutputFile(heatmap_pdf),
                mgd.OutputFile(heatmap_filt_pdf),
                mgd.OutputFile(metrics_pdf),
                mgd.OutputFile(kernel_density_pdf),
                cellids,
                config,
                args,
                params,
                params_tag,
                results_dir
            ),
            kwargs={'alignment_metrics': args['alignment_metrics']}
        )

        results = {
            'hmmcopy_metrics': helpers.format_file_yaml(hmmcopy_data),
            'segments_plot':helpers.format_file_yaml(segs_pdf),
            'bias_plot': helpers.format_file_yaml(bias_pdf),
            'filtered_heatmap_plot': helpers.format_file_yaml(heatmap_filt_pdf),
            'heatmap_plot': helpers.format_file_yaml(heatmap_pdf),
            'kde_plot': helpers.format_file_yaml(kernel_density_pdf),
            'metrics_plot': helpers.format_file_yaml(metrics_pdf)
        }

        input_datasets = {k: helpers.format_file_yaml(v) for k,v in bam_files.iteritems()}

        metadata = {
            'hmmcopy':{
                'reads_table': '/hmmcopy/reads/0',
                'parameters_table': '/hmmcopy/params/0',
                'segments_table': '/hmmcopy/segments/0',
                'metrics_table': '/hmmcopy/metrics/0',
                'hmmcopy_params_tag': params_tag,
                'hmmcopy_params': params,
                'chromosomes': config['chromosomes'],
                'ref_genome': config['ref_genome'],
                'cell_filters': config["good_cells"],
                'version': single_cell.__version__,
                'results': results,
                'containers': config['containers'],
                'input_datasets': input_datasets,
                'output_datasets': None
            }
        }

        workflow.transform(
            name='generate_meta_yaml',
            ctx=dict(mem=config['memory']['med'],
                     pool_id=config['pools']['standard'],
                     **ctx),
            func="single_cell.utils.helpers.write_to_yaml",
            args=(
                mgd.OutputFile(info_file),
                metadata
            )
        )


    return workflow
