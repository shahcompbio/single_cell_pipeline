import pypeliner
import pypeliner.managed as mgd
from single_cell.workflows.snv_allele_counts.dtypes import dtypes


def create_snv_allele_counts_for_vcf_targets_workflow(
        bam_files,
        vcf_file,
        out_file,
        memory_cfg,
        count_duplicates=False,
        min_bqual=0,
        min_mqual=0,
        vcf_to_bam_chrom_map=None,
):
    ctx = {
        'mem': memory_cfg['low'], 'num_retry': 3, 'mem_retry_increment': 2, 'ncpus': 1,
        'disk_retry_increment': 50,
    }
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'library_id', 'cell_id'),
        value=list(bam_files.keys()),
    )

    workflow.transform(
        name='get_snv_allele_counts_for_vcf_targets',
        axes=('sample_id', 'library_id', 'cell_id'),
        func="biowrappers.components.variant_calling.snv_allele_counts.tasks.get_snv_allele_counts_for_vcf_targets",
        args=(
            mgd.InputFile('tumour.bam', 'sample_id', 'library_id', 'cell_id', fnames=bam_files, extensions=['.bai']),
            mgd.InputFile(vcf_file),
            mgd.TempOutputFile('counts.csv.gz', 'sample_id', 'library_id', 'cell_id', extensions=['.yaml']),
        ),
        kwargs={
            'count_duplicates': count_duplicates,
            'min_bqual': min_bqual,
            'min_mqual': min_mqual,
            'vcf_to_bam_chrom_map': vcf_to_bam_chrom_map,
            'cell_id': mgd.Instance('cell_id'),
            'sample_id': mgd.Instance('sample_id'),
            'library_id': mgd.Instance('library_id'),
            'report_zero_count_positions': False,
            'dtypes': dtypes()['snv_allele_counts']
        }
    )

    workflow.transform(
        name='merge_snv_allele_counts',
        ctx={'mem': memory_cfg['high'], 'disk': 20},
        func="single_cell.utils.csvutils.concatenate_csv",
        args=(
            mgd.TempInputFile('counts.csv.gz', 'sample_id', 'library_id', 'cell_id'),
            mgd.OutputFile(out_file, extensions=['.yaml']),
        ),
    )

    return workflow
