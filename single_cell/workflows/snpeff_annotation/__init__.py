import pypeliner

import pypeliner.managed as mgd



def create_snpeff_annotation_workflow(
        in_vcf_file,
        out_csv_file,
        db,
        data_dir,
        split_size=int(1e3)
):
    workflow = pypeliner.workflow.Workflow(
        ctx={'num_retry': 3, 'mem_retry_increment': 2}
    )

    workflow.transform(
        name='split_vcf',
        func='single_cell.utils.vcfutils.split_vcf',
        args=(
            mgd.InputFile(in_vcf_file),
            mgd.TempOutputFile('split.vcf', 'split')
        ),
        kwargs={'lines_per_file': split_size}
    )

    workflow.transform(
        name='run_snpeff',
        axes=('split',),
        func='single_cell.workflows.snpeff_annotation.tasks.run_snpeff',
        args=(
            db,
            data_dir,
            mgd.TempInputFile('split.vcf', 'split'),
            mgd.TempOutputFile('snpeff.vcf', 'split')
        ),
        kwargs={
            'classic_mode': True
        }
    )

    workflow.transform(
        name='convert_vcf_to_csv',
        axes=('split',),
        func='single_cell.workflows.snpeff_annotation.tasks.convert_vcf_to_table',
        args=(
            mgd.TempInputFile('snpeff.vcf', 'split'),
            mgd.TempOutputFile('snpeff.csv.gz', 'split', extensions=['.yaml']),
        )
    )

    workflow.transform(
        name='concatenate_tables',
        func='single_cell.utils.csvutils.concatenate_csv',
        args=(
            mgd.TempInputFile('snpeff.csv.gz', 'split', extensions=['.yaml']),
            mgd.OutputFile(out_csv_file, extensions=['.yaml'])
        )
    )

    return workflow
