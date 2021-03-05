import pypeliner
import pypeliner.managed as mgd

from single_cell.workflows.snv_annotate.dtypes import dtypes


def create_snv_annotate_workflow(
        config,
        museq_vcf,
        strelka_vcf,
        cosmic_csv,
        dbsnp_csv,
        mappability_csv,
        snpeff_csv,
        trinuc_csv,
        memory_config,
):

    ctx = {
        'mem': memory_config['low'], 'num_retry': 3, 'mem_retry_increment': 2, 'ncpus': 1,
        'disk_retry_increment': 50,
    }
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.transform(
        name='merge_snvs',
        func='biowrappers.components.io.vcf.tasks.merge_vcfs',
        ctx=ctx,
        args=(
            [
                mgd.InputFile(museq_vcf, extensions=['.tbi', '.csi']),
                mgd.InputFile(strelka_vcf, extensions=['.tbi', '.csi']),
            ],
            mgd.TempOutputFile('all.snv.vcf')
        ),
    )

    workflow.transform(
        name='finalise_snvs',
        func="biowrappers.components.io.vcf.tasks.finalise_vcf",
        ctx=ctx,
        args=(
            mgd.TempInputFile('all.snv.vcf'),
            mgd.TempOutputFile('all.snv.vcf.gz', extensions=['.tbi', '.csi'])
        ),
    )

    workflow.subworkflow(
        name='annotate_snvs',
        axes=(),
        ctx=ctx,
        func="biowrappers.pipelines.snv_call_and_annotate.create_annotation_workflow",
        args=(
            config,
            mgd.TempInputFile('all.snv.vcf.gz', extensions=['.tbi', '.csi']),
            mgd.TempOutputFile('cosmic.csv.gz'),
            mgd.TempOutputFile('dbsnp.csv.gz'),
            mgd.TempOutputFile('mappability.csv.gz'),
            mgd.TempOutputFile('snpeff.csv.gz'),
            mgd.TempOutputFile('trinuc.csv.gz'),
        ),
        kwargs={
            'variant_type': 'snv',
        }
    )

    workflow.transform(
        name='prep_cosmic_csv',
        func='single_cell.utils.csvutils.rewrite_csv_file',
        args=(
            mgd.TempInputFile('cosmic.csv.gz'),
            mgd.OutputFile(cosmic_csv, extensions=['.yaml'])
        ),
        kwargs={'dtypes': dtypes()['snv_annotate']}
    )

    workflow.transform(
        name='prep_dbsnp_csv',
        func='single_cell.utils.csvutils.rewrite_csv_file',
        args=(
            mgd.TempInputFile('dbsnp.csv.gz'),
            mgd.OutputFile(dbsnp_csv, extensions=['.yaml'])
        ),
        kwargs={'dtypes': dtypes()['snv_annotate']}
    )

    workflow.transform(
        name='prep_mappability_csv',
        func='single_cell.utils.csvutils.rewrite_csv_file',
        args=(
            mgd.TempInputFile('mappability.csv.gz'),
            mgd.OutputFile(mappability_csv, extensions=['.yaml'])
        ),
        kwargs={'dtypes': dtypes()['snv_annotate']}
    )

    workflow.transform(
        name='prep_snpeff_csv',
        func='single_cell.utils.csvutils.rewrite_csv_file',
        args=(
            mgd.TempInputFile('snpeff.csv.gz'),
            mgd.OutputFile(snpeff_csv, extensions=['.yaml'])
        ),
        kwargs={'dtypes': dtypes()['snv_annotate']}
    )

    workflow.transform(
        name='prep_trinuc_csv',
        func='single_cell.utils.csvutils.rewrite_csv_file',
        args=(
            mgd.TempInputFile('trinuc.csv.gz'),
            mgd.OutputFile(trinuc_csv, extensions=['.yaml'])
        ),
        kwargs={'dtypes': dtypes()['snv_annotate']}
    )

    return workflow
