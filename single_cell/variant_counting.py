import os
import sys

import pypeliner.managed as mgd
from single_cell.utils import inpututils

import pypeliner


def create_snv_allele_counts_for_vcf_targets_workflow(
        bam_files,
        vcf_file,
        out_file,
        memory_cfg,
        count_duplicates=False,
        min_bqual=0,
        min_mqual=0,
        table_name='snv_allele_counts',
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
            mgd.TempOutputFile('counts.h5', 'sample_id', 'library_id', 'cell_id'),
            table_name,
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
        }
    )

    workflow.transform(
        name='merge_snv_allele_counts',
        ctx={'mem': memory_cfg['high'], 'disk': 20},
        func="biowrappers.components.io.hdf5.tasks.concatenate_tables",
        args=(
            mgd.TempInputFile('counts.h5', 'sample_id', 'library_id', 'cell_id'),
            mgd.TempOutputFile('merged_counts.h5'),
        ),
        kwargs={
            'in_memory': False,
        },
    )

    workflow.transform(
        name='convert_h5_to_csv',
        func='single_cell.utils.hdfutils.convert_hdf_to_csv',
        args=(
            mgd.TempInputFile('merged_counts.h5'),
            {
                '/snv_allele_counts': mgd.OutputFile(out_file, extensions=['.yaml']),
            }
        )
    )

    return workflow


def create_variant_counting_workflow(args):
    """ Count variant reads for multiple sets of variants across cells.
    """

    strelka_vcf, museq_vcf, tumour_cell_bams = inpututils.load_variant_counting_input(
        args['input_yaml']
    )

    counts_output = os.path.join(args['out_dir'], "counts.csv.gz")

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    config = inpututils.load_config(args)
    config = config['variant_calling']

    workflow = pypeliner.workflow.Workflow(ctx={'docker_image': config['docker']['single_cell_pipeline']})

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'library_id', 'cell_id'),
        value=list(tumour_cell_bams.keys()),
    )

    workflow.transform(
        name='merge_snvs_museq',
        func='single_cell.utils.vcfutils.merge_vcf',
        args=(
            [
                mgd.InputFile('museq.vcf', 'sample_id', 'library_id', fnames=museq_vcf,
                              extensions=['.tbi', '.csi'], axes_origin=[]),
                mgd.InputFile('strelka.vcf', 'sample_id', 'library_id', fnames=strelka_vcf,
                              extensions=['.tbi', '.csi'], axes_origin=[]),
            ],
            mgd.TempOutputFile('all.snv.vcf.gz', extensions=['.tbi', '.csi']),
            mgd.TempSpace("merge_vcf_temp")
        ),
        kwargs={'docker_image': config['docker']['vcftools']}
    )

    workflow.subworkflow(
        name='count_alleles',
        func=create_snv_allele_counts_for_vcf_targets_workflow,
        args=(
            mgd.InputFile('tumour_cells.bam', 'sample_id', 'library_id', 'cell_id', extensions=['.bai'],
                          fnames=tumour_cell_bams, axes_origin=[]),
            mgd.TempInputFile('all.snv.vcf.gz', extensions=['.tbi', '.csi']),
            mgd.OutputFile(counts_output),
            config['memory'],
        ),
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='single_cell.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args['out_dir'],
            [counts_output],
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': inpututils.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
        }
    )

    return workflow


def variant_counting_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = create_variant_counting_workflow(args)

    pyp.run(workflow)
