import os
import sys

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import inpututils


def get_file_paths(root_dir):
    return {"cohort_maf": os.path.join(root_dir, "cohort_oncogenic_filtered.maf"),
            "cohort_oncoplot": os.path.join(root_dir, "cohort_oncoplot.png"),
            "cna_table": os.path.join(root_dir, "cna_table.tsv.gz"),
            "segments": os.path.join(root_dir, "segments.tsv.gz")}


def cohort_qc_pipeline(args):
    config = inpututils.load_config(args)
    config = config["cohort_qc"]

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config['docker']['single_cell_pipeline']}
    )

    out_dir = args["out_dir"]
    api_key = args["API_key"]

    cohort, mafs, hmmcopy = inpututils.load_cohort_qc_inputs(args["input_yaml"])

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    germline_mafs = {label: data["germline_maf"] for label, data in mafs.items()}
    somatic_mafs = {label: data["somatic_maf"] for label, data in mafs.items()}
    hmmcopy_files = {label: data["hmmcopy"] for label, data in hmmcopy.items()}

    file_paths = get_file_paths(os.path.join(out_dir, cohort))

    workflow.setobj(
        obj=mgd.OutputChunks('sample_label', 'library_label'),
        value=list(hmmcopy_files.keys()),
    )

    workflow.subworkflow(
        name="classifycopynumber",
        func="single_cell.workflows.cohort_qc.cna_annotation_workflow",
        args=(
            config,
            mgd.InputFile('hmmcopy_dict', 'sample_label', 'library_label', fnames=hmmcopy_files, axes_origin=[]),
            mgd.OutputFile(file_paths["cna_table"]),
            mgd.TempOutputFile("cna_maftools_table"),
            mgd.OutputFile(file_paths["segments"]),
            config["gtf"],
        ),
    )

    workflow.subworkflow(
        name="maf_annotation_workflow",
        func="single_cell.workflows.cohort_qc.preprocess_mafs_workflow",
        args=(
            config,
            mgd.InputFile('germline_mafs_dict', 'sample_label', fnames=germline_mafs, axes_origin=[]),
            mgd.InputFile('somatic_mafs_dict', 'sample_label', fnames=somatic_mafs, axes_origin=[]),
            mgd.OutputFile(file_paths["cohort_maf"]),
            api_key
        ),
    )

    workflow.subworkflow(
        name="make_plots_and_report",
        func="single_cell.workflows.cohort_qc.create_cohort_oncoplot",
        args=(
            config,
            cohort,
            mgd.InputFile(file_paths["cohort_maf"]),
            mgd.TempInputFile("cna_maftools_table"),
            mgd.OutputFile(file_paths["cohort_oncoplot"])
        ),
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='single_cell.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args['out_dir'],
            list(file_paths.values()),
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': inpututils.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'cohort_qc'}
        }
    )
    pyp.run(workflow)
