import pypeliner
import pypeliner.managed as mgd


def create_db_annotation_workflow(
        in_vcf_file,
        out_csv_file,
        db_vcf_file,
        split_size=1e4
):
    workflow = pypeliner.workflow.Workflow(ctx=dict(mem=2, num_retry=3, mem_retry_increment=2))

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
        name='annotate_db_status',
        axes=('split',),
        func='single_cell.workflows.db_annotation.tasks.annotate_db_status',
        args=(
            db_vcf_file,
            mgd.TempInputFile('split.vcf', 'split'),
            mgd.TempOutputFile('annotated.csv.gz', 'split', extensions=['.yaml'])
        )
    )

    workflow.transform(
        name='merge_tables',
        func='single_cell.utils.csvutils.concatenate_csv',
        args=(
            mgd.TempInputFile('annotated.csv.gz', 'split', extensions=['.yaml']),
            mgd.OutputFile(out_csv_file, extensions=['.yaml'])
        )
    )

    return workflow
