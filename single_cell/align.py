'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import single_cell
import pypeliner
import pypeliner.managed as mgd
from workflows import align,alignment_metrics
from single_cell.utils import helpers


def align_workflow(args):

    config = helpers.load_config(args)
    config = config['alignment']

    baseimage = config['docker']['single_cell_pipeline']

    sampleinfo = helpers.get_sample_info(args['input_yaml'])

    cellids = helpers.get_samples(args['input_yaml'])
    bam_files, _ = helpers.get_bams(args['input_yaml'])

    lib = args["library_id"]

    outdir = os.path.join(args["out_dir"], "results", "alignment")

    info_file = os.path.join(outdir, "info.yaml")

    alignment_metrics_csv = os.path.join(outdir, '{}_alignment_metrics.csv.gz'.format(lib))
    alignment_metrics_yaml = os.path.join(outdir, '{}_alignment_metrics.yaml'.format(lib))
    gc_metrics_csv = os.path.join(outdir, '{}_gc_metrics.csv.gz'.format(lib))
    gc_metrics_yaml = os.path.join(outdir, '{}_gc_metrics.yaml'.format(lib))

    plots_dir = os.path.join(outdir,  'plots')
    plot_metrics_output = os.path.join(plots_dir, '{}_plot_metrics.pdf'.format(lib))

    workflow = pypeliner.workflow.Workflow()

    if not args["metrics_only"]:
        fastq1_files, fastq2_files = helpers.get_fastqs(args['input_yaml'])
        triminfo = helpers.get_trim_info(args['input_yaml'])
        centerinfo = helpers.get_center_info(args['input_yaml'])

        workflow.setobj(
            obj=mgd.OutputChunks('cell_id', 'lane'),
            value=fastq1_files.keys(),
        )

        workflow.subworkflow(
            name='alignment_workflow',
            ctx={'docker_image': baseimage},
            func=align.create_alignment_workflow,
            args=(
                mgd.InputFile('fastq_1', 'cell_id', 'lane', fnames=fastq1_files, axes_origin=[]),
                mgd.InputFile('fastq_2', 'cell_id', 'lane', fnames=fastq2_files, axes_origin=[]),
                mgd.OutputFile('bam_markdups', 'cell_id', fnames=bam_files, axes_origin=[], extensions=['.bai']),
                config['ref_genome'],
                config,
                args,
                triminfo,
                centerinfo,
                sampleinfo,
                cellids,
            ),
        )
    else:
        workflow.setobj(
            obj=mgd.OutputChunks('cell_id'),
            value=cellids,
        )

    workflow.subworkflow(
        name='metrics_workflow',
        func=alignment_metrics.create_alignment_metrics_workflow,
        ctx={'docker_image': baseimage},
        args=(
            mgd.InputFile('bam_markdups', 'cell_id', fnames = bam_files, axes_origin=[], extensions=['.bai']),
            mgd.OutputFile(alignment_metrics_csv),
            mgd.OutputFile(alignment_metrics_yaml),
            mgd.OutputFile(gc_metrics_csv),
            mgd.OutputFile(gc_metrics_yaml),
            mgd.OutputFile(plot_metrics_output),
            config['ref_genome'],
            config,
            args,
            sampleinfo,
            cellids,
        ),
    )

    inputs = helpers.get_fastq_files(args["input_yaml"])
    outputs = {k: helpers.format_file_yaml(v) for k,v in bam_files.iteritems()}

    metadata = {
        'alignment': {
            'name': 'alignment',
            'cell_batch_realign': args["realign"],
            'metrics_table': '/alignment/metrics',
            'gc_metrics_table': '/alignment/gc_metrics',
            'aligner': config["aligner"],
            'adapter': config["adapter"],
            'adapter2': config["adapter2"],
            'picardtools_wgsmetrics_params': config['picard_wgs_params'],
            'ref_genome': config["ref_genome"],
            'version': single_cell.__version__,
            'containers': config['docker'],
            'output_datasets': outputs,
            'input_datasets': inputs,
            'results': {
                'alignment_plots': helpers.format_file_yaml(plot_metrics_output),
            },
        }
    }

    workflow.transform(
        name='generate_meta_yaml',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.utils.helpers.write_to_yaml",
        args=(
            mgd.OutputFile(info_file),
            metadata
        )
    )

    return workflow
