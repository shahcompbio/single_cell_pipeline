import os
import pypeliner
import pypeliner.managed as mgd

from single_cell.utils import concatenate_csv

import tasks

def create_aneufinder_workflow(bam_file,
                               sample_info,
                               sample_ids,
                               config,
                               aneufinder_output,
                               aneufinder_segs_filename,
                               aneufinder_reads_filename,
                               library_id):

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )

    workflow.transform(
        name='run_aneufinder_on_individual_cells',
        ctx={'mem': config['low_mem']},
        func=tasks.run_aneufinder,
        axes=('sample_id',),
        args=(
            mgd.InputFile('bam_file', 'sample_id', fnames=bam_file),
            mgd.TempSpace('working_dir', 'sample_id', fnames=bam_file),
            mgd.InputInstance('sample_id'),
            aneufinder_output,
            mgd.TempOutputFile('segments.csv', 'sample_id'),
            mgd.TempOutputFile('reads.csv', 'sample_id'),
            mgd.TempOutputFile('dnacopy.pdf', 'sample_id'),
        ),
    )

    workflow.transform(
        name='merge_segments',
        ctx={'mem': config['low_mem']},
        func=concatenate_csv,
        args=(
            mgd.TempInputFile('segments.csv', 'sample_id'),
            mgd.OutputFile(aneufinder_segs_filename)
        ),
    )

    workflow.transform(
        name='merge_reads',
        ctx={'mem': config['low_mem']},
        func=concatenate_csv,
        args=(
            mgd.TempInputFile('reads.csv', 'sample_id'),
            mgd.OutputFile(aneufinder_reads_filename)
        )
    )

    dnacopy_pdf_output = os.path.join(aneufinder_output, 'plots', '{}_reads.pdf'.format(library_id))
    workflow.transform(
        name='merge_aneufinder_pdfs',
        ctx={'mem': config['med_mem']},
        func=tasks.merge_pdf,
        args=(
            [mgd.TempInputFile('dnacopy.pdf', 'sample_id')],
            [mgd.OutputFile(dnacopy_pdf_output)],
        )
    )

    return workflow
