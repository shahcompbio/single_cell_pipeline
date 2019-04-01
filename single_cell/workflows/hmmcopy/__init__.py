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
            mgd.TempOutputFile('reads.csv.gz', 'cell_id', 'multiplier', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile('segs.csv.gz', 'cell_id', 'multiplier', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile('params.csv.gz', 'cell_id', 'multiplier', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile('hmm_metrics.csv.gz', 'cell_id', 'multiplier', axes_origin=[], extensions=['.yaml']),
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
        func="single_cell.workflows.hmmcopy.tasks.concatenate_csv",
        args=(
            mgd.TempInputFile('reads.csv.gz', 'cell_id', 'multiplier', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile("reads.csv.gz", 'multiplier', axes_origin=[], extensions=['.yaml']),
            multipliers,
            cell_ids
        ),
        kwargs={'low_memory': True}
    )

    workflow.transform(
        name='merge_segs',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.concatenate_csv",
        args=(
            mgd.TempInputFile('segs.csv.gz', 'cell_id', 'multiplier', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile("segs.csv.gz", 'multiplier', axes_origin=[], extensions=['.yaml']),
            multipliers,
            cell_ids
        ),
        kwargs={'low_memory': True}
    )

    workflow.transform(
        name='merge_metrics',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.concatenate_csv",
        args=(
            mgd.TempInputFile('hmm_metrics.csv.gz', 'cell_id', 'multiplier', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile("hmm_metrics.csv.gz", 'multiplier', axes_origin=[], extensions=['.yaml']),
            multipliers,
            cell_ids
        ),
    )


    workflow.transform(
        name='merge_params',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.concatenate_csv",
        args=(
            mgd.TempInputFile('params.csv.gz', 'cell_id', 'multiplier', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile("params.csv.gz", 'multiplier', axes_origin=[], extensions=['.yaml']),
            multipliers,
            cell_ids
        ),
    )

    annotation_input = mgd.TempInputFile("hmm_metrics.csv.gz", "multiplier", extensions=['.yaml'])
    if alignment_metrics:
        annotation_input = mgd.TempInputFile("hmmcopy_quality_metrics.csv.gz", "multiplier", extensions=['.yaml'])
        workflow.transform(
            name="add_quality",
            ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
            axes=('multiplier',),
            func="single_cell.workflows.hmmcopy.tasks.add_quality",
            args=(
                mgd.TempInputFile("hmm_metrics.csv.gz", 'multiplier', extensions=['.yaml']),
                mgd.InputFile(alignment_metrics, extensions=['.yaml']),
                multipliers,
                mgd.TempOutputFile("hmmcopy_quality_metrics.csv.gz", "multiplier", extensions=['.yaml']),
                hmmparams['classifier_training_data'],
            ),
        )

    workflow.transform(
        name='annotate_metrics_with_info_and_clustering',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        axes=('multiplier',),
        func="single_cell.workflows.hmmcopy.tasks.annotate_metrics",
        args=(
            mgd.TempInputFile("reads.csv.gz", 'multiplier', extensions=['.yaml']),
            annotation_input,
            mgd.TempOutputFile("annotated_metrics.csv.gz", 'multiplier', extensions=['.yaml']),
            sample_info,
            cell_ids,
        ),
        kwargs={
            'chromosomes': hmmparams["chromosomes"],
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
            mgd.TempInputFile("annotated_metrics.csv.gz", 'multiplier', extensions=['.yaml']),
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
            mgd.TempInputFile("segs.csv.gz", 'multiplier', extensions=['.yaml']),
            mgd.TempInputFile("annotated_metrics.csv.gz", 'multiplier', extensions=['.yaml']),
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
            mgd.TempInputFile("annotated_metrics.csv.gz", 'multiplier', extensions=['.yaml']),
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
            mgd.TempInputFile("annotated_metrics.csv.gz", 'multiplier', extensions=['.yaml']),
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
            mgd.TempInputFile("reads.csv.gz", 'multiplier', extensions=['.yaml']),
            mgd.TempInputFile("annotated_metrics.csv.gz", 'multiplier', extensions=['.yaml']),
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
            mgd.TempInputFile("reads.csv.gz", 'multiplier', extensions=['.yaml']),
            mgd.TempInputFile("annotated_metrics.csv.gz", 'multiplier', extensions=['.yaml']),
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

    workflow.transform(
        name='finalize_reads',
        axes=('multiplier',),
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("reads.csv.gz", 'multiplier', extensions=['.yaml']),
            mgd.OutputFile("reads_final.csv.gz", 'multiplier', fnames=reads, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_segments',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        axes=('multiplier',),
        func="single_cell.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("segs.csv.gz", 'multiplier', extensions=['.yaml']),
            mgd.OutputFile("segs_final.csv.gz", 'multiplier', fnames=segs, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_metrics',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        axes=('multiplier',),
        func="single_cell.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("annotated_metrics.csv.gz", 'multiplier', extensions=['.yaml']),
            mgd.OutputFile("annotated_metrics_final.csv.gz", 'multiplier', fnames=metrics, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_params',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        axes=('multiplier',),
        func="single_cell.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("params.csv.gz", 'multiplier', extensions=['.yaml']),
            mgd.OutputFile("params_final.csv.gz", 'multiplier', fnames=params, extensions=['.yaml']),
        ),
    )

    return workflow
