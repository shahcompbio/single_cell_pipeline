'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import pypeliner.managed as mgd
from workflows import hmmcopy
from utils import helpers
import single_cell
import copy

def hmmcopy_workflow(args):
    workflow = hmmcopy_workflow(args)

    config = helpers.load_config(args)
    config = config['hmmcopy']

    cellids = helpers.get_samples(args['input_yaml'])
    bam_files, _ = helpers.get_bams(args['input_yaml'])
    lib = args['library_id']

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cellids,
    )

    for params_tag, params in config.iteritems():
        params_tag = "hmmcopy_" + params_tag

        multipliers = copy.deepcopy(params["multipliers"])
        multipliers = [0] + multipliers
        workflow.setobj(
            obj=mgd.OutputChunks('multiplier'),
            value=multipliers,
        )

        results_dir = os.path.join(args['out_dir'], 'results', params_tag)

        info_file = os.path.join(results_dir, "info.yaml")

        reads_csvs = os.path.join(
            results_dir, 'multiplier_{multiplier}', '{}_reads.csv.gz'.format(lib))
        reads_yaml = os.path.join(
            results_dir, 'multiplier_{multiplier}', '{}_reads.csv.gz.yaml'.format(lib))
        segs_csvs = os.path.join(
            results_dir, 'multiplier_{multiplier}', '{}_segments.csv.gz'.format(lib))
        segs_yaml = os.path.join(
            results_dir, 'multiplier_{multiplier}', '{}_segments.csv.gz.yaml'.format(lib))
        params_csvs = os.path.join(
            results_dir, 'multiplier_{multiplier}', '{}_params.csv.gz'.format(lib))
        params_yaml = os.path.join(
            results_dir, 'multiplier_{multiplier}', '{}_params.csv.gz.yaml'.format(lib))
        metrics_csvs = os.path.join(
            results_dir, 'multiplier_{multiplier}', '{}_metrics.csv.gz'.format(lib))
        metrics_yaml = os.path.join(
            results_dir, 'multiplier_{multiplier}', '{}_metrics.csv.gz.yaml'.format(lib))
        igv_csvs = os.path.join(
            results_dir, 'multiplier_{multiplier}', '{}_igv_segments.seg'.format(lib))


        plots_dir = os.path.join(results_dir, "plots")
        segs_pdf = os.path.join(
            plots_dir, "segments", '{}_segs.tar.gz'.format(lib))
        bias_pdf = os.path.join(plots_dir, "bias", '{}_bias.tar.gz'.format(lib))

        heatmap_filt_pdf = os.path.join(
            plots_dir, '{}_heatmap_by_ec_filtered.pdf'.format(lib))
        heatmap_pdf = os.path.join(
            plots_dir, '{}_heatmap_by_ec.pdf'.format(lib))
        metrics_pdf = os.path.join(
            plots_dir, '{}_metrics.pdf'.format(lib))
        kernel_density_pdf = os.path.join(
            plots_dir, '{}_kernel_density.pdf'.format(lib))

        baseimage = params['docker']['single_cell_pipeline']

        sample_info = helpers.get_sample_info(args["input_yaml"])

        workflow.subworkflow(
            name='hmmcopy_workflow_' + params_tag,
            func=hmmcopy.create_hmmcopy_workflow,
            ctx={'docker_image': baseimage},
            args=(
                mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_files, extensions=['.bai']),
                mgd.OutputFile("reads.csv.gz", 'multiplier', template=reads_csvs, axes_origin=[]),
                mgd.OutputFile("reads.yaml", 'multiplier', template=reads_yaml, axes_origin=[]),
                mgd.OutputFile("segs.csv.gz", 'multiplier', template=segs_csvs, axes_origin=[]),
                mgd.OutputFile("segs.yaml", 'multiplier', template=segs_yaml, axes_origin=[]),
                mgd.OutputFile("metrics.csv.gz", 'multiplier', template=metrics_csvs, axes_origin=[]),
                mgd.OutputFile("metrics.yaml", 'multiplier', template=metrics_yaml, axes_origin=[]),
                mgd.OutputFile("params.csv.gz", 'multiplier', template=params_csvs, axes_origin=[]),
                mgd.OutputFile("params.yaml", 'multiplier', template=params_yaml, axes_origin=[]),
                mgd.OutputFile("igv.seg", 'multiplier', template=igv_csvs, axes_origin=[]),
                mgd.OutputFile(segs_pdf),
                mgd.OutputFile(bias_pdf),
                mgd.OutputFile(heatmap_pdf),
                mgd.OutputFile(heatmap_filt_pdf),
                mgd.OutputFile(metrics_pdf),
                mgd.OutputFile(kernel_density_pdf),
                cellids,
                params,
                sample_info
            ),
            kwargs={'alignment_metrics': args['alignment_metrics']}
        )

        results = {
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
                'version': single_cell.__version__,
                'containers': params['docker'],
                'results': results,
                'input_datasets': input_datasets,
                'output_datasets': None
            }
        }

        workflow.transform(
            name='generate_meta_yaml',
            ctx={'mem': params['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
            func="single_cell.utils.helpers.write_to_yaml",
            args=(
                mgd.OutputFile(info_file),
                metadata
            )
        )


    return workflow
