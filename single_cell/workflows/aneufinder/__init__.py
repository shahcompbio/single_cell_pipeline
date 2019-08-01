import pypeliner.managed as mgd

import pypeliner


def create_aneufinder_workflow(bam_file,
                               cell_ids,
                               config,
                               aneufinder_results_filename,
                               aneufinder_pdf_filename,
                               ):
    baseimage = config['docker']['single_cell_pipeline']

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cell_ids,
    )

    workflow.transform(
        name='run_aneufinder_on_individual_cells',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.aneufinder.tasks.run_aneufinder",
        axes=('cell_id',),
        args=(
            mgd.InputFile('bam_file', 'cell_id', fnames=bam_file),
            mgd.TempSpace('working_dir', 'cell_id', fnames=bam_file),
            mgd.InputInstance('cell_id'),
            mgd.TempOutputFile('segments.csv', 'cell_id'),
            mgd.TempOutputFile('reads.csv', 'cell_id'),
            mgd.TempOutputFile('dnacopy.pdf', 'cell_id'),
        ),
        kwargs={'docker_image': config['docker']['aneufinder']}
    )

    workflow.transform(
        name='merge_outputs',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.aneufinder.tasks.merge_outputs_to_hdf",
        args=(
            mgd.TempInputFile('reads.csv', 'cell_id'),
            mgd.TempInputFile('segments.csv', 'cell_id'),
            mgd.OutputFile(aneufinder_results_filename),
            mgd.TempSpace("aneufinder_merge"),
        )
    )

    workflow.transform(
        name='merge_aneufinder_pdfs',
        ctx={'mem': config['memory']['med'], 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.workflows.aneufinder.tasks.merge_pdf",
        args=(
            [mgd.TempInputFile('dnacopy.pdf', 'cell_id')],
            [mgd.OutputFile(aneufinder_pdf_filename)],
        )
    )

    return workflow
