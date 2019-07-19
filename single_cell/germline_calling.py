'''
Created on Feb 22, 2018

@author: dgrewal
'''

import os
import pypeliner
import pypeliner.managed as mgd
from workflows import germline
from single_cell.utils import helpers
import single_cell

def germline_calling_workflow(args):

    config = helpers.load_config(args)
    config = config['germline_calling']

    baseimage = config['docker']['single_cell_pipeline']

    basedocker = {'docker_image': config['docker']['single_cell_pipeline']}
    vcftoolsdocker = {'docker_image': config['docker']['vcftools']}
    samtoolsdocker = {'docker_image': config['docker']['samtools']}
    snpeffdocker = {'docker_image': config['docker']['snpeff']}

    pyp = pypeliner.app.Pypeline(config=args)

    ctx = {'mem_retry_increment': 2, 'ncpus': 1, 'mem': config["memory"]['low'],
           'disk_retry_increment': 50, 'docker_image': baseimage},
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    data = helpers.load_pseudowgs_input(args['input_yaml'])
    normal_bams = data['normal_bams']
    tumour_cells = data['tumour_cells']

    if not isinstance(normal_bams, dict):
        workflow.transform(
            name="get_regions",
            func="single_cell.utils.pysamutils.get_regions_from_reference",
            ret=pypeliner.managed.OutputChunks('region'),
            args=(
                config["ref_genome"],
                config["split_size"],
                config["chromosomes"],
            )
        )
        assert '{region}' in normal_bams, 'only supports a list of files or a template on regions'
        workflow.set_filenames('normal_split.bam', 'region', template=normal_bams)
    else:
        workflow.setobj(
            obj=mgd.OutputChunks('region'),
            value=list(normal_bams.keys()),
        )
        workflow.set_filenames('normal_split.bam', 'normal_split', fnames=normal_bams)

    varcalls_dir = os.path.join(
        args['out_dir'], 'results', 'germline_calling')

    samtools_germline_vcf = os.path.join(varcalls_dir, 'raw', 'samtools_germline.vcf.gz')
    snpeff_vcf_filename = os.path.join(varcalls_dir, 'snpeff.vcf')
    normal_genotype_filename = os.path.join(varcalls_dir, 'raw', 'normal_genotype.h5')
    mappability_filename = os.path.join(varcalls_dir, 'raw', 'mappability.h5')
    counts_template = os.path.join(varcalls_dir, 'counts', 'raw', 'counts.h5')
    germline_h5_filename = os.path.join(varcalls_dir, 'germline.h5')

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=list(tumour_cells.keys()),
    )

    workflow.subworkflow(
        name='samtools_germline',
        func=germline.create_samtools_germline_workflow,
        args=(
            mgd.InputFile("normal_split.bam", "region", extensions=['.bai']),
            config['ref_genome'],
            mgd.OutputFile(samtools_germline_vcf, extensions=['.tbi']),
            config,
        ),
        kwargs={'vcftools_docker': vcftoolsdocker,
                'samtools_docker': samtoolsdocker,}
    )

    workflow.subworkflow(
        name='annotate_mappability',
        func="biowrappers.components.variant_calling.mappability.create_vcf_mappability_annotation_workflow",
        args=(
            config['databases']['mappability']['local_path'],
            mgd.InputFile(samtools_germline_vcf, extensions=['.tbi']),
            mgd.OutputFile(mappability_filename),
        ),
        kwargs={'base_docker': basedocker,
                'chromosomes': config['chromosomes']}
    )

    workflow.transform(
        name='annotate_genotype',
        func="single_cell.workflows.germline.tasks.annotate_normal_genotype",
        args=(
            mgd.InputFile(samtools_germline_vcf, extensions=['.tbi']),
            mgd.OutputFile(normal_genotype_filename),
            config["chromosomes"],
        ),
    )

    workflow.subworkflow(
        name='snpeff',
        func="biowrappers.components.variant_calling.snpeff.create_snpeff_annotation_workflow",
        args=(
            config['databases']['snpeff']['db'],
            mgd.InputFile(samtools_germline_vcf, extensions=['.tbi']),
            mgd.OutputFile(snpeff_vcf_filename),
        ),
        kwargs={
            'hdf5_output': False,
            'base_docker': basedocker,
            'vcftools_docker': vcftoolsdocker,
            'snpeff_docker': snpeffdocker,
        }
    )

    workflow.subworkflow(
        name='read_counts',
        func="single_cell.variant_calling.create_snv_allele_counts_for_vcf_targets_workflow",
        args=(
            mgd.InputFile('tumour.bam', 'cell_id', fnames=tumour_cells, extensions=['.bai']),
            mgd.InputFile(samtools_germline_vcf, extensions=['.tbi']),
            mgd.OutputFile(counts_template),
            config['memory'],
        ),
        kwargs={
            'table_name': '/germline_allele_counts',
        },
    )

    workflow.transform(
        name='build_results_file',
        func="biowrappers.components.io.hdf5.tasks.concatenate_tables",
        args=([
                mgd.InputFile(counts_template),
                mgd.InputFile(mappability_filename),
                mgd.InputFile(normal_genotype_filename),
            ],
            pypeliner.managed.OutputFile(germline_h5_filename),
        ),
        kwargs={
            'drop_duplicates': True,
        }
    )

    pyp.run(workflow)


def germline_calling_pipeline(args):

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = germline_calling_workflow(args)

    pyp.run(workflow)