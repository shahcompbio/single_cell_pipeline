import pypeliner
import pypeliner.managed as mgd
from single_cell.workflows.infer_haps.dtypes import dtypes


def infer_haps(
        bam_file,
        haplotypes_filename,
        config,
        from_tumour=False,
):
    baseimage = {'docker_image': config['docker']['single_cell_pipeline']}

    remixt_image = config['docker']['remixt']

    remixt_config = config.get('extract_seqdata', {})
    remixt_ref_data_dir = config['ref_data_dir']

    chromosomes = config['chromosomes']
    remixt_config['chromosomes'] = chromosomes

    ctx = dict(mem_retry_increment=2, disk_retry_increment=50, ncpus=1, **baseimage)
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    if isinstance(bam_file, dict):
        workflow.setobj(
            obj=mgd.OutputChunks('cell_id'),
            value=list(bam_file.keys()),
        )

        # dont parallelize over chromosomes for per cell bams
        workflow.subworkflow(
            name="extract_seqdata",
            axes=('cell_id',),
            func='remixt.workflow.create_extract_seqdata_workflow',
            ctx={'docker_image': remixt_image},
            args=(
                mgd.InputFile(
                    'bam_markdups', 'cell_id', fnames=bam_file, extensions=['.bai']
                ),
                mgd.TempOutputFile('seqdata_cell.h5', 'cell_id'),
                remixt_config,
                remixt_ref_data_dir,
            ),
            kwargs={'no_parallelism': True}
        )
        workflow.transform(
            name='merge_all_seqdata',
            func="remixt.seqdataio.merge_overlapping_seqdata",
            ctx={'docker_image': remixt_image},
            args=(
                mgd.TempOutputFile('seqdata_file.h5'),
                mgd.TempInputFile("seqdata_cell.h5", "cell_id"),
                config["chromosomes"]
            ),
        )
    else:
        workflow.subworkflow(
            name='extract_seqdata',
            func='remixt.workflow.create_extract_seqdata_workflow',
            ctx={'disk': 150, 'docker_image': remixt_image},
            args=(
                mgd.InputFile(bam_file, extensions=['.bai']),
                mgd.TempOutputFile('seqdata_file.h5'),
                remixt_config,
                remixt_ref_data_dir,
            ),
        )

    workflow.setobj(
        obj=mgd.OutputChunks('chromosome'),
        value=chromosomes,
    )

    if from_tumour:
        func = 'remixt.analysis.haplotype.infer_snp_genotype_from_tumour'
    else:
        func = 'remixt.analysis.haplotype.infer_snp_genotype_from_normal'

    workflow.transform(
        name='infer_snp_genotype',
        axes=('chromosome',),
        ctx={'mem': 16, 'docker_image': remixt_image},
        func=func,
        args=(
            mgd.TempOutputFile('snp_genotype.tsv', 'chromosome'),
            mgd.TempInputFile('seqdata_file.h5'),
            mgd.InputInstance('chromosome'),
            config,
        ),
    )

    workflow.transform(
        name='infer_haps',
        axes=('chromosome',),
        ctx={'mem': 16, 'docker_image': remixt_image},
        func='remixt.analysis.haplotype.infer_haps',
        args=(
            mgd.TempOutputFile('haplotypes.tsv', 'chromosome'),
            mgd.TempInputFile('snp_genotype.tsv', 'chromosome'),
            mgd.InputInstance('chromosome'),
            mgd.TempSpace('haplotyping', 'chromosome'),
            remixt_config,
            remixt_ref_data_dir,
        ),
    )

    workflow.transform(
        name='merge_haps',
        ctx={'mem': 16, 'docker_image': remixt_image},
        func='remixt.utils.merge_tables',
        args=(
            mgd.TempOutputFile('haplotypes_merged.tsv'),
            mgd.TempInputFile('haplotypes.tsv', 'chromosome'),
        )
    )

    workflow.transform(
        name='finalize_csv',
        ctx={'mem': 16},
        func='single_cell.utils.csvutils.rewrite_csv_file',
        args=(
            mgd.TempInputFile('haplotypes_merged.tsv'),
            mgd.OutputFile(haplotypes_filename, extensions=['.yaml']),
        ),
        kwargs={
            'write_header': True,
            'dtypes': dtypes()['haplotypes']
        },
    )

    return workflow
