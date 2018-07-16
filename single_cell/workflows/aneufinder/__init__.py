import os
import pypeliner
import pypeliner.managed as mgd


def create_aneufinder_workflow(bam_file,
                               cell_ids,
                               config,
                               aneufinder_output,
                               aneufinder_results_filename,
                               library_id):

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cell_ids,
    )

    workflow.transform(
        name='run_aneufinder_on_individual_cells',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'mem_retry_increment':2,
            'ncpus': 1},
        func="single_cell.workflows.aneufinder.tasks.run_aneufinder",
        axes=('cell_id',),
        args=(
            mgd.InputFile('bam_file', 'cell_id', fnames=bam_file),
            mgd.TempSpace('working_dir', 'cell_id', fnames=bam_file),
            mgd.InputInstance('cell_id'),
            aneufinder_output,
            mgd.TempOutputFile('segments.csv', 'cell_id'),
            mgd.TempOutputFile('reads.csv', 'cell_id'),
            mgd.TempOutputFile('dnacopy.pdf', 'cell_id'),
        ),
    )

    workflow.transform(
        name='merge_outputs',
        ctx={
            'mem': config["memory"]['low'],
            'pool_id': config['pools']['standard'],
            'mem_retry_increment':2,
            'ncpus': 1},
        func="single_cell.workflows.aneufinder.tasks.merge_outputs_to_hdf",
        args=(
            mgd.TempInputFile('reads.csv', 'cell_id'),
            mgd.TempInputFile('segments.csv', 'cell_id'),
            mgd.OutputFile(aneufinder_results_filename),
            mgd.TempSpace("aneufinder_merge"),
        )
    )

    dnacopy_pdf_output = os.path.join(
        aneufinder_output,
        'plots',
        '{}_reads.pdf'.format(library_id))
    workflow.transform(
        name='merge_aneufinder_pdfs',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'mem_retry_increment':2,
            'ncpus': 1},
        func="single_cell.workflows.aneufinder.tasks.merge_pdf",
        args=(
            [mgd.TempInputFile('dnacopy.pdf', 'cell_id')],
            [mgd.OutputFile(dnacopy_pdf_output)],
        )
    )

    return workflow
