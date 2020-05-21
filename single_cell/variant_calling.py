'''
Created on Feb 22, 2018

@author: dgrewal
'''

import os
import sys

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import inpututils
from single_cell.workflows import mutationseq
from single_cell.workflows import snv_annotate
from single_cell.workflows import strelka


def get_file_paths(root_dir):
    data = {
        'museq_vcf': os.path.join(root_dir, 'museq.vcf.gz'),
        'cosmic_csv': os.path.join(root_dir, 'snv_cosmic_status.csv.gz'),
        'dbsnp_csv': os.path.join(root_dir, 'snv_dbsnp_status.csv.gz'),
        'mappability_csv': os.path.join(root_dir, 'snv_mappability.csv.gz'),
        'snpeff_csv': os.path.join(root_dir, 'snv_snpeff.csv.gz'),
        'museq_csv': os.path.join(root_dir, 'snv_museq.csv.gz'),
        'strelka_csv': os.path.join(root_dir, 'snv_strelka.csv.gz'),
        'trinuc_csv': os.path.join(root_dir, 'snv_trinuc.csv.gz'),
        'strelka_indel': os.path.join(root_dir, 'strelka_indel.vcf.gz'),
        'strelka_snv': os.path.join(root_dir, 'strelka_snv.vcf.gz'),
    }

    return data


def variant_calling_workflow(args):
    config = inpututils.load_config(args)
    config = config['variant_calling']

    normal_bams, tumour_bams = inpututils.load_variant_calling_input(args['input_yaml'])

    filepaths = get_file_paths(args['out_dir'])

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    ctx = {
        'ncpus': 1,
        'mem_retry_increment': 2,
        'disk_retry_increment': 50,
        'mem': config["memory"]['low'],
        'docker_image': config['docker']['single_cell_pipeline']
    }
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('region'),
        value=list(normal_bams.keys()),
    )
    workflow.subworkflow(
        name='museq',
        func=mutationseq.create_museq_workflow,
        args=(
            mgd.InputFile('normal_regions.bam', 'region', extensions=['.bai'], fnames=normal_bams),
            mgd.InputFile('tumour_regions.bam', 'region', extensions=['.bai'], fnames=tumour_bams),
            mgd.OutputFile(filepaths['museq_vcf'], extensions=['.tbi', '.csi']),
            mgd.OutputFile(filepaths['museq_csv'], extensions=['.tbi', '.csi']),
            config,
        ),
    )

    workflow.subworkflow(
        name='strelka',
        func=strelka.create_strelka_workflow,
        args=(
            mgd.InputFile('normal_regions.bam', 'region', extensions=['.bai'], fnames=normal_bams),
            mgd.InputFile('tumour_regions.bam', 'region', extensions=['.bai'], fnames=tumour_bams),
            config['ref_genome'],
            mgd.OutputFile(filepaths['strelka_indel'], extensions=['.tbi', '.csi']),
            mgd.OutputFile(filepaths['strelka_snv'], extensions=['.tbi', '.csi']),
            mgd.OutputFile(filepaths['strelka_csv'], extensions=['.yaml']),
            config,
        ),
        kwargs={
            "chromosomes": config["chromosomes"],
            "use_depth_thresholds": config['use_depth_thresholds']
        }
    )

    workflow.subworkflow(
        name='annotate_snv_vcf_files',
        func=snv_annotate.create_snv_annotate_workflow,
        args=(
            config,
            mgd.InputFile(filepaths['museq_vcf'], extensions=['.tbi', '.csi']),
            mgd.InputFile(filepaths['strelka_snv'], extensions=['.tbi', '.csi']),
            mgd.OutputFile(filepaths['cosmic_csv'], extensions=['.yaml']),
            mgd.OutputFile(filepaths['dbsnp_csv'], extensions=['.yaml']),
            mgd.OutputFile(filepaths['mappability_csv'], extensions=['.yaml']),
            mgd.OutputFile(filepaths['snpeff_csv'], extensions=['.yaml']),
            mgd.OutputFile(filepaths['trinuc_csv'], extensions=['.yaml']),
            config['docker'],
            config['memory']
        )
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='single_cell.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args['out_dir'],
            list(filepaths.values()),
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': inpututils.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'variant_calling'}
        }
    )

    return workflow


def variant_calling_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = variant_calling_workflow(args)

    pyp.run(workflow)
