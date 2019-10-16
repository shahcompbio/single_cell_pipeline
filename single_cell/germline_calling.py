'''
Created on Feb 22, 2018

@author: dgrewal
'''

import os
import sys

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import inpututils
from single_cell.workflows import germline


def get_output_files(outdir):
    data = {
        'samtools_germline_vcf': os.path.join(outdir, 'samtools_germline.vcf.gz'),
        'snpeff_vcf_filename': os.path.join(outdir, 'snpeff.vcf.gz'),
        'normal_genotype_filename': os.path.join(outdir, 'normal_genotype..vcf.gz'),
        'mappability_filename': os.path.join(outdir, 'mappability.csv.gz'),

    }

    return data


def germline_calling_workflow(args):
    config = inpututils.load_config(args)
    config = config['germline_calling']

    vcftoolsdocker = {'docker_image': config['docker']['vcftools']}
    samtoolsdocker = {'docker_image': config['docker']['samtools']}
    snpeffdocker = {'docker_image': config['docker']['snpeff']}

    normal_bams = inpututils.load_germline_data(args['input_yaml'])

    varcalls_meta = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')
    out_files = get_output_files(args['out_dir'])

    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config['docker']['single_cell_pipeline']}
    )

    workflow.setobj(
        obj=mgd.OutputChunks('region'),
        value=list(normal_bams.keys()),
    )

    workflow.subworkflow(
        name='samtools_germline',
        func=germline.create_samtools_germline_workflow,
        args=(
            mgd.InputFile("normal_split.bam", "region", extensions=['.bai'], fnames=normal_bams),
            config['ref_genome'],
            mgd.OutputFile(out_files['samtools_germline_vcf'], extensions=['.tbi']),
            config,
        ),
        kwargs={'vcftools_docker': vcftoolsdocker,
                'samtools_docker': samtoolsdocker, }
    )

    workflow.subworkflow(
        name='annotate_mappability',
        func="biowrappers.components.variant_calling.mappability.create_vcf_mappability_annotation_workflow",
        args=(
            config['databases']['mappability']['local_path'],
            mgd.InputFile(out_files['samtools_germline_vcf'], extensions=['.tbi']),
            mgd.OutputFile(out_files['mappability_filename']),
        ),
        kwargs={'chromosomes': config['chromosomes']}
    )

    workflow.transform(
        name='annotate_genotype',
        func="single_cell.workflows.germline.tasks.annotate_normal_genotype",
        args=(
            mgd.InputFile(out_files['samtools_germline_vcf'], extensions=['.tbi']),
            mgd.OutputFile(out_files['normal_genotype_filename']),
            config["chromosomes"],
        ),
    )

    workflow.subworkflow(
        name='snpeff',
        func="biowrappers.components.variant_calling.snpeff.create_snpeff_annotation_workflow",
        args=(
            config['databases']['snpeff']['db'],
            mgd.InputFile(out_files['samtools_germline_vcf'], extensions=['.tbi']),
            mgd.OutputFile(out_files['snpeff_vcf_filename']),
        ),
        kwargs={
            'hdf5_output': False,
            'vcftools_docker': vcftoolsdocker,
            'snpeff_docker': snpeffdocker,
        }
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='single_cell.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args['out_dir'],
            list(out_files.values()),
            mgd.OutputFile(varcalls_meta)
        ),
        kwargs={
            'input_yaml_data': inpututils.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'germline_calling'}
        }
    )

    return workflow


def germline_calling_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = germline_calling_workflow(args)

    pyp.run(workflow)
