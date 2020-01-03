import logging
import os

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import inpututils


def create_sv_genotyper_workflow(
        tumour_cell_bams,
        lumpy_calls,
        destruct_calls,
        out_dir,
        config
):
    annotations = [
        "AO", "AP", "AS", "ASC", "DP", "GQ", "QA",
        "QR", "RO", "RP", "RS", "SQ", "GL", "AB"
    ]
    output_template = {ann: os.path.join(out_dir, '{sample_id}', '{library_id}', ann + ".csv.gz") for ann in
                       annotations}

    ctx = {
        'mem': 8, 'num_retry': 3, 'mem_retry_increment': 2, 'ncpus': 1,
        'disk_retry_increment': 50,
    }

    sv_genotyping_config = config['sv_genotyping']
    reference_genome = sv_genotyping_config['ref_genome']

    error_str = '''
  ___ __  __ ___  ___  ___  ___  __  __  ___  _  _  _____  _    _      ___  ___    _  _____  _   _  ___  ___
 | __|\ \/ /| _ \| __|| _ \|_ _||  \/  || __|| \| ||_   _|/_\  | |    | __|| __|  /_\|_   _|| | | || _ \| __|
 | _|  >  < |  _/| _| |   / | | | |\/| || _| | .` |  | | / _ \ | |__  | _| | _|  / _ \ | |  | |_| ||   /| _|
 |___|/_/\_\|_|  |___||_|_\|___||_|  |_||___||_|\_|  |_|/_/ \_\|____| |_|  |___|/_/ \_\|_|   \___/ |_|_\|___|
 '''
    logging.getLogger("SV_GENOTYPING").warning(error_str)

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'library_id', 'cell_id'),
        value=list(tumour_cell_bams.keys()),
    )

    workflow.subworkflow(
        name='genotype',
        func='single_cell.workflows.sv_genotyping.create_sv_genotyping_workflow',
        axes=('sample_id', 'library_id',),
        ctx={'docker_image': sv_genotyping_config['docker']['single_cell_pipeline']},
        args=(
            reference_genome,
            mgd.InputFile('bam_input', 'sample_id', 'library_id', 'cell_id', extensions=['.bai'],
                          axes_origin=[], fnames=tumour_cell_bams),
            mgd.InputFile('lumpy_csv_gz', 'sample_id', 'library_id', axes_origin=[], fnames=lumpy_calls),
            mgd.InputFile('destruct_csv_gz', 'sample_id', 'library_id', axes_origin=[],
                          fnames=destruct_calls),
            mgd.OutputFile('AO', 'sample_id', 'library_id', template=output_template["AO"]),
            mgd.OutputFile('AP', 'sample_id', 'library_id', template=output_template["AP"]),
            mgd.OutputFile('AS', 'sample_id', 'library_id', template=output_template["AS"]),
            mgd.OutputFile('ASC', 'sample_id', 'library_id', template=output_template["ASC"]),
            mgd.OutputFile('DP', 'sample_id', 'library_id', template=output_template["DP"]),
            mgd.OutputFile('GQ', 'sample_id', 'library_id', template=output_template["GQ"]),
            mgd.OutputFile('QA', 'sample_id', 'library_id', template=output_template["QA"]),
            mgd.OutputFile('QR', 'sample_id', 'library_id', template=output_template["QR"]),
            mgd.OutputFile('RO', 'sample_id', 'library_id', template=output_template["RO"]),
            mgd.OutputFile('RP', 'sample_id', 'library_id', template=output_template["RP"]),
            mgd.OutputFile('RS', 'sample_id', 'library_id', template=output_template["RS"]),
            mgd.OutputFile('SQ', 'sample_id', 'library_id', template=output_template["SQ"]),
            mgd.OutputFile('GL', 'sample_id', 'library_id', template=output_template["GL"]),
            mgd.OutputFile('AB', 'sample_id', 'library_id', template=output_template["AB"]),
            sv_genotyping_config
        )
    )

    return workflow


def sv_genotyping_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    lumpy_sv, destruct_sv, tumour_bams = inpututils.load_sv_genotyper_input(args['input_yaml'])

    out_dir = args['out_dir']

    config = inpututils.load_config(args)

    workflow = create_sv_genotyper_workflow(
        tumour_bams,
        lumpy_sv,
        destruct_sv,
        out_dir,
        config
    )

    pyp.run(workflow)
