import pypeliner
import pypeliner.managed as mgd


def create_snv_annotate_workflow(
        config,
        museq_vcf,
        strelka_vcf,
        mappability_csv,
        snpeff_csv,
        trinuc_csv,
        additional_csv,
        memory_config,
):
    ctx = {
        'mem': memory_config['low'], 'num_retry': 3, 'mem_retry_increment': 2, 'ncpus': 1,
        'disk_retry_increment': 50,
    }
    split_size = config['split_size']

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
        name='snpeff_annotation',
        func="single_cell.workflows.snpeff_annotation.create_snpeff_annotation_workflow",
        args=(
            mgd.TempInputFile('all.snv.vcf.gz', extensions=['.tbi', '.csi']),
            mgd.OutputFile(snpeff_csv, extensions=['.yaml']),
            config['databases']['snpeff']['db'],
            config['databases']['snpeff']['path'],
        )
    )

    workflow.subworkflow(
        name='trinuc_annotation',
        func="single_cell.workflows.trinuc_annotation.create_trinuc_annotation_workflow",
        args=(
            mgd.TempInputFile('all.snv.vcf.gz', extensions=['.tbi', '.csi']),
            mgd.OutputFile(trinuc_csv, extensions=['.yaml']),
            config['ref_genome'],
        ),
        kwargs={'split_size': split_size}
    )

    workflow.subworkflow(
        name='mappability_annotation',
        func="single_cell.workflows.mappability_annotation.create_mappability_annotation_workflow",
        args=(
            mgd.TempInputFile('all.snv.vcf.gz', extensions=['.tbi', '.csi']),
            mgd.OutputFile(mappability_csv, extensions=['.yaml']),
            config['databases']['mappability']['path'],
        ),
        kwargs={'split_size': split_size}
    )

    for k, v in config['databases']['additional_databases'].items():
        workflow.subworkflow(
            name='{}_status'.format(k),
            func='single_cell.workflows.db_annotation.create_db_annotation_workflow',
            ctx=dict(mem=4, mem_retry_increment=2),
            args=(
                mgd.TempInputFile('all.snv.vcf.gz', extensions=['.tbi', '.csi']),
                mgd.OutputFile(additional_csv[k], extensions=['.yaml']),
                v['path'],
            ),
            kwargs={'split_size': split_size}
        )

    return workflow
