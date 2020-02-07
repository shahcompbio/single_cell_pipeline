'''
Created on Jul 6, 2017

@author: dgrewal
'''
import pypeliner.managed as mgd

import pypeliner


def create_hmmcopy_workflow(
        bam_file, reads, segs, metrics, params, igv_seg_filename,
        segs_pdf, bias_pdf, plot_heatmap_ec_output,
        plot_metrics_output,
        plot_kernel_density_output, hmmcopy_data_tar,
        cell_ids, hmmparams, sample_info
):
    chromosomes = hmmparams["chromosomes"]

    baseimage = hmmparams['docker']['single_cell_pipeline']
    hmmcopy_docker = hmmparams['docker']['hmmcopy']

    ctx = {'mem': 7, 'ncpus': 1, 'docker_image': baseimage, 'mem_retry_factor': 1}

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

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
            mgd.TempOutputFile('reads.csv.gz', 'cell_id', extensions=['.yaml']),
            mgd.TempOutputFile('segs.csv.gz', 'cell_id', extensions=['.yaml']),
            mgd.TempOutputFile('params.csv.gz', 'cell_id', extensions=['.yaml']),
            mgd.TempOutputFile('hmm_metrics.csv.gz', 'cell_id', extensions=['.yaml']),
            mgd.TempOutputFile('hmm_data.tar.gz', 'cell_id'),
            mgd.InputInstance('cell_id'),
            hmmparams,
            mgd.TempSpace('hmmcopy_temp', 'cell_id'),
            hmmcopy_docker
        ),
    )

    workflow.transform(
        name='merge_reads',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.concatenate_csv",
        args=(
            mgd.TempInputFile('reads.csv.gz', 'cell_id', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile('reads_merged.csv.gz', extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='add_mappability_bool',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.get_mappability_col",
        args=(
            mgd.TempInputFile('reads_merged.csv.gz', extensions=['.yaml']),
            mgd.OutputFile(reads, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='merge_segs',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.concatenate_csv",
        args=(
            mgd.TempInputFile('segs.csv.gz', 'cell_id', axes_origin=[], extensions=['.yaml']),
            mgd.OutputFile(segs, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='merge_metrics',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.concatenate_csv",
        args=(
            mgd.TempInputFile('hmm_metrics.csv.gz', 'cell_id', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile("hmm_metrics.csv.gz", extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='merge_params',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.concatenate_csv",
        args=(
            mgd.TempInputFile('params.csv.gz', 'cell_id', axes_origin=[], extensions=['.yaml']),
            mgd.OutputFile(params, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='get_max_cn',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.get_max_cn",
        ret=mgd.TempOutputObj('max_cn'),
        args=(
            mgd.InputFile(reads, extensions=['.yaml']),
        )
    )

    workflow.transform(
        name='hmmcopy_plots',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.plot_hmmcopy",
        axes=('cell_id',),
        args=(
            mgd.TempInputFile('reads.csv.gz', 'cell_id', axes_origin=[], extensions=['.yaml']),
            mgd.TempInputFile('segs.csv.gz', 'cell_id', axes_origin=[], extensions=['.yaml']),
            mgd.TempInputFile('params.csv.gz', 'cell_id', axes_origin=[], extensions=['.yaml']),
            mgd.TempInputFile('hmm_metrics.csv.gz', 'cell_id', axes_origin=[], extensions=['.yaml']),
            hmmparams['ref_genome'],
            mgd.TempOutputFile('segments.png', 'cell_id', axes_origin=[]),
            mgd.TempOutputFile('bias.png', 'cell_id', axes_origin=[]),
            mgd.InputInstance('cell_id'),
        ),
        kwargs={
            'num_states': hmmparams['num_states'],
            'sample_info': mgd.TempInputObj('sampleinfo', 'cell_id'),
            'max_cn': mgd.TempInputObj("max_cn")
        }
    )

    workflow.transform(
        name='annotate_metrics_with_info_and_clustering',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.add_clustering_order",
        args=(
            mgd.InputFile(reads, extensions=['.yaml']),
            mgd.TempInputFile("hmm_metrics.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(metrics, extensions=['.yaml']),
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
            mgd.InputFile(metrics, extensions=['.yaml']),
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
            mgd.InputFile(segs, extensions=['.yaml']),
            mgd.InputFile(metrics, extensions=['.yaml']),
            mgd.OutputFile(igv_seg_filename),
            hmmparams,
        )
    )

    workflow.transform(
        name='plot_metrics',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.plot_metrics",
        args=(
            mgd.InputFile(metrics, extensions=['.yaml']),
            mgd.OutputFile(plot_metrics_output),
            'QC pipeline metrics',
        )
    )

    workflow.transform(
        name='plot_kernel_density',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.hmmcopy.tasks.plot_kernel_density",
        args=(
            mgd.InputFile(metrics, extensions=['.yaml']),
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
            mgd.InputFile(reads, extensions=['.yaml']),
            mgd.InputFile(metrics, extensions=['.yaml']),
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
        name='merge_hmmcopy_data_tars',
        ctx={'mem': hmmparams['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.utils.helpers.tar_files",
        args=(
            mgd.TempInputFile('hmm_data.tar.gz', 'cell_id', axes_origin=[]),
            mgd.OutputFile(hmmcopy_data_tar),
            mgd.TempSpace("merge_tarballs")
        ),

    )

    return workflow
