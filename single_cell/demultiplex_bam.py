'''
Created on November 01, 2018

@author: pashaa
'''
import os
import glob
import single_cell
import pypeliner
import pypeliner.managed as mgd
from workflows import demultiplex_bam
from single_cell.utils import helpers



def demultiplex_bam_workflow(workflow, args):
    #fileConfig("logging_config.ini")

    config = helpers.load_config(args)

    out_dir = os.path.join(args["out_dir"], os.path.basename(args['bam'][:-4]))

    info_file = os.path.join(out_dir, "info.yaml")

    ctx = {'mem_retry_increment': 2, 'ncpus': 1}
    ctx.update(helpers.get_container_ctx(config['containers'], 'single_cell_pipeline'))

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


    inputs = {'filename': args['bam']}

    outputs = {'dirname': out_dir}

    metadata = {
        'demultiplex_bam': {
            'name': 'demultiplex_bam',
            'version': single_cell.__version__,
            'containers': config['containers'],
            'output_datasets': outputs,
            'input_datasets': inputs,
            'results': {
                'demultiplexed_fastqs': glob.glob(os.path.join(out_dir, "*.fq"))
            },
        }
    }

    workflow.transform(
        name='generate_meta_yaml',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        func="single_cell.utils.helpers.write_to_yaml",
        args=(
            mgd.OutputFile(info_file),
            metadata
        )
    )

    return workflow
