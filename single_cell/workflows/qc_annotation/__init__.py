'''
Created on Jul 6, 2017

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd


def create_qc_annotation_workflow(
        hmmcopy_metrics, alignment_metrics, gc_metrics, merged_metrics, qc_report,
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

    workflow.transform(
        name='generate_qc_report',
        func="single_cell.workflows.qc_annotation.tasks.generate_qc_report",
        args=(
            mgd.TempSpace("QC_report_singlecellpipeline"),
            config['reference_gc'],
            mgd.InputFile(merged_metrics, extensions=['.yaml']),
            mgd.InputFile(gc_metrics, extensions=['.yaml']),
            mgd.OutputFile(qc_report)
        )
    )

    return workflow
