'''
Created on Jul 6, 2017

@author: dgrewal
'''
import copy
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers


def create_hmmcopy_workflow(
        bam_file, reads, reads_yaml, segs, segs_yaml, metrics, metrics_yaml,
        params, params_yaml, igv_seg_filename, segs_pdf, bias_pdf,
        plot_heatmap_ec_output, plot_heatmap_ec_filt_output,
        plot_metrics_output, plot_kernel_density_output, cell_ids,
        hmmparams, sample_info, alignment_metrics=None
):
    chromosomes = hmmparams["chromosomes"]

    multipliers = copy.deepcopy(hmmparams["multipliers"])
    multipliers = [0] + multipliers

    igv_seg_filename = dict([(mult, igv_seg_filename[mult])
                             for mult in multipliers])
    reads = dict([(mult, reads[mult])
                  for mult in multipliers])
    reads_yaml = dict([(mult, reads_yaml[mult])
                       for mult in multipliers])
    segs = dict([(mult, segs[mult])
                 for mult in multipliers])
    segs_yaml = dict([(mult, segs_yaml[mult])
                      for mult in multipliers])
    metrics = dict([(mult, metrics[mult])
                    for mult in multipliers])
    metrics_yaml = dict([(mult, metrics_yaml[mult])
                         for mult in multipliers])
    params = dict([(mult, params[mult])
                   for mult in multipliers])
    params_yaml = dict([(mult, params_yaml[mult])
                        for mult in multipliers])

    baseimage = hmmparams['docker']['single_cell_pipeline']
    hmmcopy_docker = hmmparams['docker']['hmmcopy']

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('multiplier'),
        value=multipliers,
    )

    allfiles = [(cell, mult) for cell in cell_ids for mult in multipliers]
    workflow.setobj(
        obj=mgd.OutputChunks('cell_id', 'multiplier'),
        value=allfiles,
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('sampleinfo', 'cell_id', axes_origin=[]),
        value=sample_info)

    workflow.transform(
        name='run_hmmcopy',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.run_hmmcopy",
        axes=('cell_id',),
        args=(
            mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_file, extensions=['.bai']),
            mgd.TempOutputFile('reads.csv.gz', 'cell_id', 'multiplier', axes_origin=[]),
            mgd.TempOutputFile('segs.csv.gz', 'cell_id', 'multiplier', axes_origin=[]),
            mgd.TempOutputFile('params.csv.gz', 'cell_id', 'multiplier', axes_origin=[]),
            mgd.TempOutputFile('hmm_metrics.csv.gz', 'cell_id', 'multiplier', axes_origin=[]),
            mgd.TempOutputFile('segments.png', 'cell_id', axes_origin=[]),
            mgd.TempOutputFile('bias.png', 'cell_id', axes_origin=[]),
            mgd.InputInstance('cell_id'),
            hmmparams['ref_genome'],
            hmmparams,
            multipliers,
            mgd.TempSpace('hmmcopy_temp', 'cell_id'),
            mgd.TempInputObj('sampleinfo', 'cell_id'),
            hmmcopy_docker
        ),
    )

    workflow.transform(
        name='merge_reads',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        axes=('multiplier',),
        func="single_cell.workflows.hmmcopy.tasks.concatenate_csv_single",
        args=(
            mgd.TempInputFile('reads.csv.gz', 'cell_id', 'multiplier', axes_origin=[]),
            mgd.OutputFile("reads.csv.gz", 'multiplier', fnames=reads),
            mgd.OutputFile("reads.yaml", 'multiplier', fnames=reads_yaml),
            mgd.InputInstance('multiplier'),
            cell_ids
        ),
        kwargs={'low_memory': True, 'quick': True}
    )

    workflow.transform(
        name='merge_segs',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.concatenate_csv",
        args=(
            mgd.TempInputFile('segs.csv.gz', 'cell_id', 'multiplier', axes_origin=[]),
            mgd.OutputFile("segs.csv.gz", 'multiplier', fnames=segs, axes_origin=[]),
            mgd.OutputFile("segs.yaml", 'multiplier', fnames=segs_yaml, axes_origin=[]),
            multipliers,
            cell_ids
        ),
        kwargs={'low_memory': True, 'quick': True}
    )

    workflow.transform(
        name='merge_metrics',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.concatenate_csv",
        args=(
            mgd.TempInputFile('hmm_metrics.csv.gz', 'cell_id', 'multiplier', axes_origin=[]),
            mgd.TempOutputFile("hmm_metrics.csv.gz", 'multiplier', axes_origin=[]),
            None,
            multipliers,
            cell_ids
        ),
    )

    workflow.transform(
        name='merge_params',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.concatenate_csv",
        args=(
            mgd.TempInputFile('params.csv.gz', 'cell_id', 'multiplier', axes_origin=[]),
            mgd.OutputFile("params.csv.gz", 'multiplier', fnames=params, axes_origin=[]),
            mgd.OutputFile("params.yaml", 'multiplier', fnames=params_yaml, axes_origin=[]),
            multipliers,
            cell_ids
        ),
    )

    annotation_input = mgd.TempInputFile("hmm_metrics.csv.gz", "multiplier")
    if alignment_metrics:
        annotation_input = mgd.TempInputFile("hmmcopy_quality_metrics.csv.gz", "multiplier")
        workflow.transform(
            name="add_quality",
            ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
            axes=('multiplier',),
            func="single_cell.workflows.hmmcopy.tasks.add_quality",
            args=(
                mgd.TempInputFile("hmm_metrics.csv.gz", 'multiplier'),
                mgd.InputFile(alignment_metrics),
                multipliers,
                mgd.TempOutputFile("hmmcopy_quality_metrics.csv.gz", "multiplier"),
                hmmparams['classifier_training_data'],
            ),
        )

    workflow.transform(
        name='annotate_metrics_with_info_and_clustering',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        axes=('multiplier',),
        func="single_cell.workflows.hmmcopy.tasks.annotate_metrics",
        args=(
            mgd.InputFile("reads.csv.gz", 'multiplier', fnames=reads),
            annotation_input,
            mgd.OutputFile("annotated_metrics.csv.gz", 'multiplier', fnames=metrics),
            sample_info,
            cell_ids,
        ),
        kwargs={
            'chromosomes': hmmparams["chromosomes"],
            'yamlfile': mgd.OutputFile(
                "annotated_metrics.yaml", 'multiplier', fnames=metrics_yaml
            ),
        }
    )

    workflow.transform(
        name='merge_hmm_copy_plots',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.merge_pdf",
        args=(
            [
                mgd.TempInputFile('segments.png', 'cell_id'),
                mgd.TempInputFile('bias.png', 'cell_id'),
            ],
            [
                mgd.OutputFile(segs_pdf),
                mgd.OutputFile(bias_pdf),
            ],
            mgd.InputFile("annotated_metrics.csv.gz", 'multiplier', fnames=metrics),
            None,
            mgd.TempSpace("hmmcopy_plot_merge_temp"),
            ['segments', 'bias']
        )
    )

    workflow.transform(
        name='create_igv_seg',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        axes=('multiplier',),
        func="single_cell.workflows.hmmcopy.tasks.create_igv_seg",
        args=(
            mgd.InputFile("segs.csv.gz", 'multiplier', fnames=segs),
            mgd.InputFile("annotated_metrics.csv.gz", 'multiplier', fnames=metrics),
            mgd.OutputFile("igv.seg", 'multiplier', fnames=igv_seg_filename),
            hmmparams,
            mgd.InputInstance("multiplier")
        )
    )

    workflow.transform(
        name='plot_metrics',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        axes=('multiplier',),
        func="single_cell.workflows.hmmcopy.tasks.plot_metrics",
        args=(
            mgd.InputFile("annotated_metrics.csv.gz", 'multiplier', fnames=metrics),
            mgd.TempOutputFile("metrics.pdf", 'multiplier'),
            'QC pipeline metrics',
            mgd.InputInstance('multiplier')
        )
    )

    workflow.transform(
        name='merge_metrics_plots',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.utils.pdfutils.merge_pdfs",
        args=(
            mgd.TempInputFile("metrics.pdf", 'multiplier'),
            mgd.OutputFile(plot_metrics_output)
        )
    )

    workflow.transform(
        name='plot_kernel_density',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        axes=('multiplier',),
        func="single_cell.workflows.hmmcopy.tasks.plot_kernel_density",
        args=(
            mgd.InputFile("annotated_metrics.csv.gz", 'multiplier', fnames=metrics),
            mgd.TempOutputFile("kde.pdf", 'multiplier'),
            ',',
            'mad_neutral_state',
            'QC pipeline metrics',
            mgd.InputInstance('multiplier')
        )
    )

    workflow.transform(
        name='merge_kde_plots',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.utils.pdfutils.merge_pdfs",
        args=(
            mgd.TempInputFile("kde.pdf", 'multiplier'),
            mgd.OutputFile(plot_kernel_density_output),
        )
    )

    workflow.transform(
        name='plot_heatmap_ec',
        axes=('multiplier',),
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.plot_pcolor",
        args=(
            mgd.InputFile("reads.csv.gz", 'multiplier', fnames=reads),
            mgd.InputFile("annotated_metrics.csv.gz", 'multiplier', fnames=metrics),
            mgd.TempOutputFile('plot_heatmap_ec_output.pdf', 'multiplier'),
            mgd.InputInstance('multiplier')
        ),
        kwargs={
            'plot_title': 'QC pipeline metrics',
            'column_name': 'state',
            'plot_by_col': 'experimental_condition',
            'color_by_col': 'cell_call',
            'chromosomes': chromosomes,
            'max_cn': hmmparams['num_states'],
            'scale_by_cells': False,
            'mappability_threshold': hmmparams["map_cutoff"]
        }
    )

    workflow.transform(
        name='merge_heatmaps',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.utils.pdfutils.merge_pdfs",
        args=(
            mgd.TempInputFile('plot_heatmap_ec_output.pdf', 'multiplier'),
            mgd.OutputFile(plot_heatmap_ec_output),
        )
    )

    workflow.transform(
        name='plot_heatmap_ec_filtered',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        axes=('multiplier',),
        func="single_cell.workflows.hmmcopy.tasks.plot_pcolor",
        args=(
            mgd.InputFile("reads.csv.gz", 'multiplier', fnames=reads),
            mgd.InputFile("annotated_metrics.csv.gz", 'multiplier', fnames=metrics),
            mgd.TempOutputFile('plot_heatmap_ec_filt_output.pdf', 'multiplier'),
            mgd.InputInstance('multiplier')
        ),
        kwargs={
            'plot_title': 'QC pipeline metrics',
            'column_name': 'state',
            'plot_by_col': 'experimental_condition',
            'color_by_col': 'cell_call',
            'chromosomes': chromosomes,
            'max_cn': hmmparams['num_states'],
            'scale_by_cells': False,
            'cell_filters': hmmparams["good_cells"],
            'mappability_threshold': hmmparams["map_cutoff"]
        }
    )

    workflow.transform(
        name='merge_filtered_heatmaps',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.utils.pdfutils.merge_pdfs",
        args=(
            mgd.TempInputFile('plot_heatmap_ec_filt_output.pdf', 'multiplier'),
            mgd.OutputFile(plot_heatmap_ec_filt_output),
        )
    )

    return workflow
