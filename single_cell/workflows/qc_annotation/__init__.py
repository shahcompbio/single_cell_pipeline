'''
Created on Jul 6, 2017

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd


def create_qc_annotation_workflow(
        hmmcopy_metrics, alignment_metrics, merged_metrics,
        sample_info, config, multipliers, cell_ids
):
    ctx = {'docker_image': config['docker']['single_cell_pipeline']}

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

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
        name="add_quality",
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        axes=('multiplier',),
        func="single_cell.workflows.qc_annotation.tasks.add_quality",
        args=(
            mgd.InputFile("hmm_metrics.csv.gz", 'multiplier', extensions=['.yaml'], fnames=hmmcopy_metrics),
            mgd.InputFile(alignment_metrics, extensions=['.yaml']),
            multipliers,
            mgd.TempOutputFile("hmmcopy_quality_metrics.csv.gz", "multiplier", extensions=['.yaml']),
            config['classifier_training_data'],
            mgd.TempSpace("hmmcopy_classify_tempdir", 'multiplier')
        ),
    )

    workflow.transform(
        name='annotate_metrics_with_info_and_clustering',
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        axes=('multiplier',),
        func="single_cell.workflows.qc_annotation.tasks.annotate_metrics",
        args=(
            mgd.TempInputFile("hmmcopy_quality_metrics.csv.gz", "multiplier", extensions=['.yaml']),
            mgd.TempOutputFile("annotated_metrics.csv.gz", 'multiplier', extensions=['.yaml']),
            sample_info,
            cell_ids,
        ),
    )

    workflow.transform(
        name='merge_alignment_hmmcopy_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        axes=('multiplier',),
        func="single_cell.workflows.qc_annotation.tasks.merge_metrics",
        args=(
            mgd.TempOutputFile("annotated_metrics.csv.gz", 'multiplier', extensions=['.yaml']),
            mgd.InputFile(alignment_metrics, extenstions=['.yaml']),
            mgd.OutputFile('all_metrics.csv.gz', 'multiplier', extensions=['.yaml'], fnames=merged_metrics)
        )
    )

    return workflow
