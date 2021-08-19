import pypeliner
import pypeliner.managed as mgd


def create_mappability_annotation_workflow(
        in_vcf_file,
        out_csv_file,
        mappability_file,
        split_size=1e4
):
    workflow = pypeliner.workflow.Workflow(
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2}
    )

    workflow.transform(
        name="get_regions",
        func="single_cell.workflows.mappability_annotation.tasks.get_vcf_regions",
        ret=mgd.TempOutputObj('regions_obj', 'regions'),
        args=(
            mgd.InputFile(in_vcf_file, extensions=['.tbi']),
            int(split_size),
        ),
    )

    workflow.transform(
        name='annotate_db_status',
        axes=('regions',),
        func='single_cell.workflows.mappability_annotation.tasks.get_mappability',
        args=(
            mappability_file,
            mgd.InputFile(in_vcf_file, extensions=['.tbi']),
            mgd.TempOutputFile('mappability.csv.gz', 'regions', extensions=['.yaml'])
        ),
        kwargs={
            'region': mgd.TempInputObj('regions_obj', 'regions'),
        },
    )

    workflow.transform(
        name='merge_tables',
        func='single_cell.utils.csvutils.concatenate_csv',
        args=(
            mgd.TempInputFile('mappability.csv.gz', 'regions', extensions=['.yaml']),
            mgd.OutputFile(out_csv_file, extensions=['.yaml'])
        )
    )

    return workflow
