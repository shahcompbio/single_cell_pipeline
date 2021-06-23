import pypeliner
import pypeliner.managed as mgd


def create_trinuc_annotation_workflow(
        in_vcf_file,
        out_csv_file,
        ref_genome,
        split_size=int(1e4),
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
        name='annotate_db_status',
        axes=('split',),
        func='single_cell.workflows.trinuc_annotation.tasks.get_tri_nucelotide_context',
        args=(
            ref_genome,
            mgd.TempInputFile('split.vcf', 'split'),
            mgd.TempOutputFile('tri_nucleotide_context.csv.gz', 'split', extensions=['.yaml']),
        )
    )

    workflow.transform(
        name='merge_tables',
        func='single_cell.utils.csvutils.concatenate_csv',
        args=(
            mgd.TempInputFile('tri_nucleotide_context.csv.gz', 'split', extensions=['.yaml']),
            mgd.OutputFile(out_csv_file, extensions=['.yaml'])
        )
    )

    return workflow
