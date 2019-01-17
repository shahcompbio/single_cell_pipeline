'''
Created on Jul 6, 2017

@author: dgrewal
'''
import copy
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers

def create_hmmcopy_workflow(
        bam_file, hmmcopy_data, igv_seg_filename, segs_pdf, bias_pdf,
        plot_heatmap_ec_output,
        plot_heatmap_ec_filt_output, plot_metrics_output,
        plot_kernel_density_output, cell_ids, args,
        hmmparams, sample_info, alignment_metrics=None):

    chromosomes = hmmparams["chromosomes"]

    multipliers = copy.deepcopy(hmmparams["multipliers"])
    multipliers = [0] + multipliers

    baseimage = hmmparams['docker']['single_cell_pipeline']
    hmmcopy_docker = hmmparams['docker']['hmmcopy']

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cell_ids,
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
            mgd.TempOutputFile('reads.h5', 'cell_id'),
            mgd.TempOutputFile('segs.h5', 'cell_id'),
            mgd.TempOutputFile('params.h5', 'cell_id'),
            mgd.TempOutputFile('hmm_metrics.h5', 'cell_id'),
            mgd.TempOutputFile('segments.png', 'cell_id'),
            mgd.TempOutputFile('bias.png', 'cell_id'),
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
        func="single_cell.workflows.hmmcopy.tasks.merge_hdf_files_on_disk",
        args=(
            mgd.TempInputFile('reads.h5', 'cell_id'),
            mgd.TempOutputFile("reads.h5"),
            multipliers,
            'hmmcopy/reads'
        ),
        kwargs={
            'dtypes': {'valid': bool, 'ideal': bool, 'state':float, 'multiplier':float}
        }
    )

    workflow.transform(
        name='merge_segs',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.merge_hdf_files_on_disk",
        args=(
            mgd.TempInputFile('segs.h5', 'cell_id'),
            mgd.TempOutputFile("segments.h5"),
            multipliers,
            'hmmcopy/segments'
        ),
        kwargs={
            'dtypes': {'end': float, 'median': float}
        }
    )

    workflow.transform(
        name='merge_metrics',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.merge_hdf_files_in_memory",
        args=(
            mgd.TempInputFile('hmm_metrics.h5', 'cell_id'),
            mgd.TempOutputFile("hmmcopy_metrics.h5"),
            multipliers,
            'hmmcopy/metrics'
        ),
        kwargs={
            'dtypes': {'mad_neutral_state': float}
        }
    )

    workflow.transform(
        name='merge_params',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.merge_hdf_files_in_memory",
        args=(
            mgd.TempInputFile('params.h5', 'cell_id'),
            mgd.TempOutputFile("params.h5"),
            multipliers,
            'hmmcopy/params'
        ),
    )

    annotation_input = 'hmmcopy_metrics.h5'
    if alignment_metrics:
        annotation_input = 'hmmcopy_quality_metrics.h5'
        workflow.transform(
            name="add_quality",
            ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
            func="single_cell.workflows.hmmcopy.tasks.add_quality",
            args=(
                mgd.TempInputFile('hmmcopy_metrics.h5'),
                mgd.InputFile(alignment_metrics),
                multipliers,
                mgd.TempOutputFile("hmmcopy_quality_metrics.h5"),
                hmmparams['classifier_training_data'],
            ),
        )

    workflow.transform(
        name='annotate_metrics_with_info_and_clustering',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.annotate_metrics",
        args=(
            mgd.TempInputFile('reads.h5'),
            mgd.TempInputFile(annotation_input),
            mgd.TempOutputFile("annotated_metrics.h5"),
            sample_info,
            cell_ids,
            multipliers,
        ),
        kwargs={'chromosomes': hmmparams["chromosomes"]}
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
            mgd.TempInputFile("annotated_metrics.h5"),
            None,
            mgd.TempSpace("hmmcopy_plot_merge_temp"),
            ['segments','bias']
        )
    )

    workflow.transform(
        name='create_igv_seg',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.create_igv_seg",
        args=(
            mgd.TempInputFile("segments.h5"),
            mgd.TempInputFile("annotated_metrics.h5"),
            mgd.OutputFile(igv_seg_filename),
            hmmparams
        )
    )

    workflow.transform(
        name='plot_metrics',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.plot_metrics",
        args=(
            mgd.TempInputFile("annotated_metrics.h5"),
            mgd.OutputFile(plot_metrics_output),
            mgd.TempSpace("plot_metrics_temp"),
            'QC pipeline metrics',
            multipliers
        )
    )

    workflow.transform(
        name='plot_kernel_density',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.plot_kernel_density",
        args=(
            mgd.TempInputFile('annotated_metrics.h5'),
            mgd.OutputFile(plot_kernel_density_output),
            mgd.TempSpace("hmmcopy_kde_plot_temp"),
            ',',
            'mad_neutral_state',
            'QC pipeline metrics',
            multipliers
        )
    )

    workflow.transform(
        name='plot_heatmap_ec',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.plot_pcolor",
        args=(
            mgd.TempInputFile('reads.h5'),
            mgd.TempInputFile('annotated_metrics.h5'),
            mgd.OutputFile(plot_heatmap_ec_output),
            mgd.TempSpace("heatmap_ec_temp"),
            multipliers
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
            mgd.TempInputFile('reads.h5'),
            mgd.TempInputFile('annotated_metrics.h5'),
            mgd.OutputFile(plot_heatmap_ec_filt_output),
            mgd.TempSpace("heatmap_ec_filt_temp"),
            multipliers
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
        name='merge_all_hdf5_stores',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.merge_tables",
        args=(
            mgd.TempInputFile("reads.h5"),
            mgd.TempInputFile("segments.h5"),
            mgd.TempInputFile("annotated_metrics.h5"),
            mgd.TempInputFile("params.h5"),
            mgd.TempOutputFile("hmmcopy_precast.h5"),
            cell_ids
        )
    )

    workflow.transform(
        name='cast_h5',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.utils.hdfutils.cast_h5_file",
        args=(
            mgd.TempInputFile("hmmcopy_precast.h5"),
            mgd.OutputFile(hmmcopy_data),
        )
    )

    return workflow
