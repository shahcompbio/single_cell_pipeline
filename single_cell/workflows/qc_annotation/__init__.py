'''
Created on Jul 6, 2017

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd


def create_qc_annotation_workflow(
        hmmcopy_metrics, alignment_metrics, merged_metrics,
        config
):
    ctx = {'docker_image': config['docker']['single_cell_pipeline']}

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.transform(
        name="add_quality",
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        func="single_cell.workflows.qc_annotation.tasks.add_quality",
        args=(
            mgd.InputFile(hmmcopy_metrics, extensions=['.yaml']),
            mgd.InputFile(alignment_metrics, extensions=['.yaml']),
            mgd.TempOutputFile("hmmcopy_quality_metrics.csv.gz", extensions=['.yaml']),
            config['classifier_training_data'],
            mgd.TempSpace("hmmcopy_classify_tempdir")
        ),
    )

    workflow.transform(
        name='merge_alignment_hmmcopy_metrics',
        ctx={'mem': config['memory']['med'], 'ncpus': 1},
        func="single_cell.workflows.qc_annotation.tasks.merge_metrics",
        args=(
            mgd.TempInputFile("hmmcopy_quality_metrics.csv.gz", extensions=['.yaml']),
            mgd.InputFile(alignment_metrics, extenstions=['.yaml']),
            mgd.OutputFile(merged_metrics, extensions=['.yaml'])
        )
    )

    return workflow
