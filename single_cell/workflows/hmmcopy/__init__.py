'''
Created on Jul 6, 2017

@author: dgrewal
'''
import copy
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers


def create_hmmcopy_workflow(
        bam_file, reads, segs, metrics, params, igv_seg_filename,
        segs_pdf, bias_pdf, plot_heatmap_ec_output,
        plot_heatmap_ec_filt_output, plot_metrics_output,
        plot_kernel_density_output, cell_ids, hmmparams, sample_info
):
    chromosomes = hmmparams["chromosomes"]

    multipliers = copy.deepcopy(hmmparams["multipliers"])
    multipliers = [0] + multipliers

    baseimage = hmmparams['docker']['single_cell_pipeline']
    hmmcopy_docker = hmmparams['docker']['hmmcopy']

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('multiplier'),
        value=multipliers,
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('sampleinfo', 'cell_id', axes_origin=[]),
        value=sample_info)

    allfiles = [(cell, mult) for cell in cell_ids for mult in multipliers]
    workflow.setobj(
        obj=mgd.OutputChunks('cell_id', 'multiplier'),
        value=allfiles,
    )

    workflow.transform(
        name='run_hmmcopy',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.run_hmmcopy",
        axes=('cell_id',),
        args=(
            mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_file, extensions=['.bai']),
            mgd.TempOutputFile('reads.csv.gz', 'cell_id', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile('segs.csv.gz', 'cell_id', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile('params.csv.gz', 'cell_id', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile('hmm_metrics.csv.gz', 'cell_id', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile('segments.png', 'cell_id', axes_origin=[]),
            mgd.TempOutputFile('bias.png', 'cell_id', axes_origin=[]),
            mgd.InputInstance('cell_id'),
            hmmparams['ref_genome'],
            hmmparams,
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
            mgd.TempInputFile('reads.csv.gz', 'cell_id', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile("reads.csv.gz", axes_origin=[], extensions=['.yaml']),
        ),
        kwargs={'low_memory': True}
    )

    workflow.transform(
        name='merge_segs',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.concatenate_csv",
        args=(
            mgd.TempInputFile('segs.csv.gz', 'cell_id', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile("segs.csv.gz", axes_origin=[], extensions=['.yaml']),
        ),
        kwargs={'low_memory': True}
    )

    workflow.transform(
        name='merge_metrics',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.concatenate_csv",
        args=(
            mgd.TempInputFile('hmm_metrics.csv.gz', 'cell_id', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile("hmm_metrics.csv.gz", axes_origin=[], extensions=['.yaml']),
        ),
    )


    workflow.transform(
        name='merge_params',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.concatenate_csv",
        args=(
            mgd.TempInputFile('params.csv.gz', 'cell_id', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile("params.csv.gz", axes_origin=[], extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='annotate_metrics_with_info_and_clustering',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.add_clustering_order",
        args=(
            mgd.TempInputFile("reads.csv.gz", extensions=['.yaml']),
            mgd.TempInputFile("hmm_metrics.csv.gz", extensions=['.yaml']),
            mgd.TempOutputFile("annotated_metrics.csv.gz", extensions=['.yaml']),
        ),
        kwargs={
            'chromosomes': hmmparams["chromosomes"],
            'sample_info': sample_info
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
            mgd.TempInputFile("annotated_metrics.csv.gz", extensions=['.yaml']),
            None,
            mgd.TempSpace("hmmcopy_plot_merge_temp"),
            ['segments', 'bias']
        )
    )

    workflow.transform(
        name='create_igv_seg',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.create_igv_seg",
        args=(
            mgd.TempInputFile("segs.csv.gz", extensions=['.yaml']),
            mgd.TempInputFile("annotated_metrics.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(igv_seg_filename),
            hmmparams,
        )
    )

    workflow.transform(
        name='plot_metrics',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.plot_metrics",
        args=(
            mgd.TempInputFile("annotated_metrics.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(plot_metrics_output),
            'QC pipeline metrics',
        )
    )

    workflow.transform(
        name='plot_kernel_density',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.plot_kernel_density",
        args=(
            mgd.TempInputFile("annotated_metrics.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(plot_kernel_density_output),
            ',',
            'mad_neutral_state',
            'QC pipeline metrics',
        )
    )


    workflow.transform(
        name='plot_heatmap_ec',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.plot_pcolor",
        args=(
            mgd.TempInputFile("reads.csv.gz", extensions=['.yaml']),
            mgd.TempInputFile("annotated_metrics.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(plot_heatmap_ec_output),
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
        name='plot_heatmap_ec_filtered',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.plot_pcolor",
        args=(
            mgd.TempInputFile("reads.csv.gz", extensions=['.yaml']),
            mgd.TempInputFile("annotated_metrics.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(plot_heatmap_ec_filt_output),
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
        name='finalize_reads',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("reads.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(reads, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_segments',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("segs.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(segs, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_metrics',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("annotated_metrics.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(metrics, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_params',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("params.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(params, extensions=['.yaml']),
        ),
    )

    return workflow
