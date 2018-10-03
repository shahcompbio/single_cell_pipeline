import os
import pypeliner
import pypeliner.managed as mgd
import single_cell
from single_cell.utils import helpers

def create_aneufinder_workflow(bam_file,
                               cell_ids,
                               config,
                               aneufinder_output,
                               aneufinder_results_filename,
                               library_id,
                               meta_yaml):

    ctx = {'mem_retry_increment': 2, 'ncpus': 1}
    docker_ctx = helpers.get_container_ctx(config['containers'], 'single_cell_pipeline')
    ctx.update(docker_ctx)

    aneufinder_docker = helpers.get_container_ctx(config['containers'], 'aneufinder', docker_only=True)

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cell_ids,
    )

    workflow.transform(
        name='run_aneufinder_on_individual_cells',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
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
        kwargs={'docker_config': aneufinder_docker}
    )

    workflow.transform(
        name='merge_outputs',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
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
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        func="single_cell.workflows.aneufinder.tasks.merge_pdf",
        args=(
            [mgd.TempInputFile('dnacopy.pdf', 'cell_id')],
            [mgd.OutputFile(dnacopy_pdf_output)],
        )
    )


    results = {
        'aneufinder_plot': helpers.format_file_yaml(dnacopy_pdf_output),
        'aneufinder_data':helpers.format_file_yaml(aneufinder_results_filename),
    }

    input_datasets = {k: helpers.format_file_yaml(v) for k,v in bam_file.iteritems()}

    metadata = {
        'aneufinder':{
            'reads_table': '/aneufinder/reads',
            'segments_table': '/aneufinder/segments/',
            'chromosomes': config['chromosomes'],
            'ref_genome': config['ref_genome'],
            'version': single_cell.__version__,
            'results': results,
            'containers': config['containers'],
            'input_datasets': input_datasets,
            'output_datasets': None
        }
    }

    workflow.transform(
        name='generate_meta_yaml',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        func="single_cell.utils.helpers.write_to_yaml",
        args=(
            mgd.OutputFile(meta_yaml),
            metadata
        )
    )


    return workflow
