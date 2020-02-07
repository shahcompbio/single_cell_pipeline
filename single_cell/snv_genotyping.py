import os
import sys

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import inpututils


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
        func='single_cell.workflows.snv_allele_counts.create_snv_allele_counts_for_vcf_targets_workflow',
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
            'metadata': {'type': 'snv_genotyping'}
        }
    )

    return workflow


def snv_genotyping_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = create_variant_counting_workflow(args)

    pyp.run(workflow)
