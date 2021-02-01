import os
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import inpututils
import sys

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

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    #inputs
    cohort, mafs, hmmcopy = inpututils.load_cohort_qc_inputs(args["input_yaml"])

    out_dir = args["out_dir"]
    tmp_dir = args["tmpdir"]
    api_key = args["API_key"]
    gtf = config["gtf"]

    germline_mafs = {label: data["germline_maf"] for label, data in mafs.items()}
    somatic_mafs = {label: data["somatic_maf"] for label, data in mafs.items()}
    hmmcopy_files = {label: data["hmmcopy"] for label, data in hmmcopy.items()}

    #outputs
    filepaths = get_file_paths( os.path.join(out_dir, cohort) )

    cna_cbioportal_table = filepaths["cna_table"]
    segments = filepaths["segments"]
    cohort_maf_oncogenic_filtered = filepaths["cohort_maf"]
    cohort_oncoplot = filepaths["cohort_oncoplot"]


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
            mgd.OutputFile(cna_cbioportal_table),
            mgd.TempOutputFile("cna_maftools_table"),
            mgd.OutputFile(segments),
            gtf,
        ),
    )

    workflow.subworkflow(
        name="maf_annotation_workflow",
        func="single_cell.workflows.cohort_qc.preprocess_mafs_workflow",
        args=(
            config,
            mgd.InputFile('germline_mafs_dict',  'sample_label', fnames=germline_mafs, axes_origin=[]),
            mgd.InputFile('somatic_mafs_dict',  'sample_label', fnames=somatic_mafs, axes_origin=[]),
            mgd.OutputFile(cohort_maf_oncogenic_filtered),
            api_key
        ),
    )

    workflow.subworkflow(
        name="make_plots_and_report",
        func="single_cell.workflows.cohort_qc.create_cohort_oncoplot",
        args=(
            config,
            cohort,
            out_dir,
            mgd.InputFile(cohort_maf_oncogenic_filtered),
            mgd.TempInputFile("cna_maftools_table"),
            mgd.OutputFile(cohort_oncoplot)
        ),
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
            'metadata': {'type': 'cohort_qc'}
        }
    )
    pyp.run(workflow)
