'''
Created on November 01, 2018

@author: pashaa
'''
import os
import single_cell
import pypeliner
import pypeliner.managed as mgd
from workflows import demultiplex_bam
from single_cell.utils import helpers
from logging.config import fileConfig



def demultiplex_bam_workflow(workflow, args):
    #fileConfig("logging_config.ini")

    config = helpers.load_config(args)

    out_dir = os.path.join(args["out_dir"], "results", args['bam'])

    info_file = os.path.join(out_dir, "info.yaml")

    # ctx = {'mem_retry_increment': 2, 'ncpus': 1}
    # ctx.update(helpers.get_container_ctx(config['containers'], 'single_cell_pipeline'))

    workflow.subworkflow(
        name='create_demultiplex_bam_workflow',
        axes=(),
        func=demultiplex_bam.create_demultiplex_bam_workflow,
        args=(
            pypeliner.managed.InputFile(args['bam']),
            config,
            out_dir,
            args['barcode_csv'],
        ),
    )




    # TODO fix input and output dict
    #inputs = helpers.get_fastq_files(args["input_yaml"])

    #outputs = {k: helpers.format_file_yaml(v) for k,v in bam_files.iteritems()}
    # inputs = {}
    # outputs = {}
    #
    # metadata = {
    #     'demultiplex_bam': {
    #         'name': 'demultiplex_bam',
    #         'version': single_cell.__version__,
    #         'containers': config['containers'],
    #         'output_datasets': outputs,
    #         'input_datasets': inputs,
    #         'results': {
    #             'demultiplexed_fastqs': helpers.format_file_yaml(out_dir)
    #         },
    #     }
    # }
    #
    # workflow.transform(
    #     name='generate_meta_yaml',
    #     ctx=dict(mem=config['memory']['med'],
    #              pool_id=config['pools']['standard'],
    #              **ctx),
    #     func="single_cell.utils.helpers.write_to_yaml",
    #     args=(
    #         mgd.OutputFile(info_file),
    #         metadata
    #     )
    # )

    return workflow
