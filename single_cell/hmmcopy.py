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
import pypeliner

def hmmcopy_workflow(args):
    workflow = pypeliner.workflow.Workflow()

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
            results_dir, '{0}_multiplier{1}_reads.csv.gz'.format(lib, '{multiplier}'))
        segs_csvs = os.path.join(
            results_dir, '{0}_multiplier{1}_segments.csv.gz'.format(lib, '{multiplier}'))
        params_csvs = os.path.join(
            results_dir, '{0}_multiplier{1}_params.csv.gz'.format(lib, '{multiplier}'))
        metrics_csvs = os.path.join(
            results_dir, '{0}_multiplier{1}_metrics.csv.gz'.format(lib, '{multiplier}'))
        igv_csvs = os.path.join(
            results_dir, '{0}_multiplier{1}_igv_segments.seg'.format(lib, '{multiplier}'))

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

        workflow.subworkflow(
            name='hmmcopy_workflow_' + params_tag,
            func=hmmcopy.create_hmmcopy_workflow,
            ctx={'docker_image': baseimage},
            args=(
                mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_files, extensions=['.bai']),
                mgd.OutputFile("reads.csv.gz", 'multiplier', template=reads_csvs, axes_origin=[]),
                mgd.OutputFile("segs.csv.gz", 'multiplier', template=segs_csvs, axes_origin=[]),
                mgd.OutputFile("metrics.csv.gz", 'multiplier', template=metrics_csvs, axes_origin=[]),
                mgd.OutputFile("params.csv.gz", 'multiplier', template=params_csvs, axes_origin=[]),
                mgd.OutputFile("igv.seg", 'multiplier', template=igv_csvs, axes_origin=[]),
                mgd.OutputFile(segs_pdf),
                mgd.OutputFile(bias_pdf),
                mgd.OutputFile(heatmap_pdf),
                mgd.OutputFile(heatmap_filt_pdf),
                mgd.OutputFile(metrics_pdf),
                mgd.OutputFile(kernel_density_pdf),
                cellids,
                params,
            ),
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
