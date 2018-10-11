'''
Created on Feb 22, 2018

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
from workflows import merge_bams
from single_cell.utils import helpers
import single_cell

def merge_bams_workflow(workflow, args):


    input_yaml = args["input_yaml"]
    output_template = args["merged_bam_template"]

    info_file = os.path.join(args["out_dir"], 'results','merge_bams', "info.yaml")
    config = helpers.load_config(args)
    bam_files, bai_files  = helpers.get_bams(input_yaml)
    cellids = helpers.get_samples(input_yaml)

    wgs_bam_template = output_template
    wgs_bai_template = wgs_bam_template + ".bai"

    ctx = {'mem_retry_increment': 2, 'ncpus': 1}
    ctx.update(helpers.get_container_ctx(config['containers'], 'single_cell_pipeline'))


    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cellids,
    )

    workflow.transform(
        name="get_regions",
        ctx=dict(mem=2, pool_id=config['pools']['standard'], **ctx),
        func="single_cell.utils.pysamutils.get_regions_from_reference",
        ret=pypeliner.managed.TempOutputObj('region'),
        args=(
              config["ref_genome"],
              config["split_size"],
              config["chromosomes"],
        )
    )

    workflow.subworkflow(
        name="wgs_merge_workflow",
        func=merge_bams.create_merge_bams_workflow,
        args=(
            mgd.InputFile('bam_markdups', 'cell_id', fnames=bam_files, extensions=['.bai']),
            mgd.OutputFile("merged_bam", "region", axes_origin=[], template=wgs_bam_template, extensions=['.bai']),
            cellids,
            config,
            mgd.TempInputObj("region"),
        )
    )

    workflow.transform(
        name="get_files",
        ctx=dict(mem=2, pool_id=config['pools']['standard'], **ctx),
        func='single_cell.utils.helpers.resolve_template',
        ret=pypeliner.managed.TempOutputObj('outputs'),
        args=(
            pypeliner.managed.TempInputObj('region'),
            wgs_bam_template,
            'region'
        )

    )

    inputs = {k: helpers.format_file_yaml(v) for k,v in bam_files.iteritems()}

    metadata = {
        'merge_bams': {
            'name': 'merge_bams',
            'ref_genome': config["ref_genome"],
            'version': single_cell.__version__,
            'containers': config['containers'],
            'output_datasets': pypeliner.managed.TempInputObj('outputs'),
            'input_datasets': inputs,
            'results': None
        }
    }

    workflow.transform(
        name='generate_meta_yaml',
        ctx=dict(mem=2, pool_id=config['pools']['standard'], **ctx),
        func="single_cell.utils.helpers.write_to_yaml",
        args=(
            mgd.OutputFile(info_file),
            metadata
        )
    )


    return workflow


