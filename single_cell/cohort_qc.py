import os
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import inpututils
import sys


def get_cbioportal_paths(root_dir):
    """Get cbioportal output paths.

    Args:
        root_dir ([str]): [path to out_dir]

    Returns:
        [dict]: [labeled output paths]
    """
    if not os.path.exists(os.path.join(root_dir, "cbioportal")):
        os.makedirs(os.path.join(root_dir, "cbioportal"))

    filtered_germline_maf = os.path.join(
        root_dir, "cbioportal", "filtered_germline.maf"
    )
    annotated_somatic_maf = os.path.join(
        root_dir, "cbioportal", "annotated_somatic.maf"
    )
    cna_table = os.path.join(
        root_dir, "cbioportal",  "cna_table.tsv"
    )
    segments = os.path.join(
        root_dir, "cbioportal",  "segments.tsv"
    )

    return {
        "filtered_germline_maf": filtered_germline_maf,
        "annotated_somatic_maf": annotated_somatic_maf,
        "cna_table": cna_table,
        "segments": segments
    }


def get_maftools_paths(root_dir):
    """Get maftools output paths.

    Args:
        root_dir ([str]): [path to out_dir]

    Returns:
        [dict]: [labeled output paths]
    """

    if not os.path.exists(os.path.join(root_dir, "maftools")):
        os.makedirs(os.path.join(root_dir, "maftools"))

    cohort_oncoplot = os.path.join(
        root_dir, "maftools", "cohort_oncoplot.maf"
    )
    maftools_maf = os.path.join(
        root_dir,  "maftools", "maftools_maf.maf"
    )
    maftools_cna = os.path.join(
        root_dir, "maftools",  "maftools_cna.tsv"
    )

    return {
        "cohort_oncoplot": cohort_oncoplot,
        "maftools_maf": maftools_maf,
        "maftools_cna": maftools_cna
    }


def cohort_qc_pipeline(args):
    """Process maf, run classify copynumber, make plots.

    Args:
        args ([dict]): [pipeline arguments]
    """
    config = inpututils.load_config(args)
    config = config["cohort_qc"]

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config['docker']['single_cell_pipeline']}
    )

    out_dir = args["out_dir"]
    api_key = args["API_key"]

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    # inputs
    mafs, hmmcopy_filenames, hmmcopy_samples = inpututils.load_cohort_qc_inputs(
        args["input_yaml"]
    )

    germline_mafs = {
        label: data["germline_maf"] for label, data in mafs.items()
    }
    somatic_mafs = {
        label: data["somatic_maf"] for label, data in mafs.items()
    }

    # outputs
    cbiofile_paths = get_cbioportal_paths(out_dir)
    maftools_filepaths = get_maftools_paths(out_dir)

    workflow.setobj(
        obj=mgd.OutputChunks('sample_label', 'library_label'),
        value=list(hmmcopy_filenames.keys()),
    )

    workflow.subworkflow(
        name="classifycopynumber",
        func="single_cell.workflows.cohort_qc.cna_annotation_workflow",
        args=(
            config,
            mgd.InputFile(
                'hmmcopy_dict', 'sample_label', 'library_label',
                fnames=hmmcopy_filenames, axes_origin=[]
            ),
            hmmcopy_samples,
            mgd.OutputFile(cbiofile_paths["cna_table"]),
            mgd.OutputFile(maftools_filepaths["maftools_cna"]),
            mgd.OutputFile(cbiofile_paths["segments"]),
            config["gtf"],
        ),
    )

    workflow.subworkflow(
        name="maf_annotation_workflow",
        func="single_cell.workflows.cohort_qc.preprocess_mafs_workflow",
        args=(
            config,
            mgd.InputFile(
                'germline_mafs_dict',  'sample_label',
                fnames=germline_mafs, axes_origin=[]
            ),
            mgd.InputFile(
                'somatic_mafs_dict',  'sample_label',
                fnames=somatic_mafs, axes_origin=[]
            ),
            mgd.OutputFile(cbiofile_paths["filtered_germline_maf"]),
            mgd.OutputFile(cbiofile_paths["annotated_somatic_maf"]),
            api_key
        ),
    )
    workflow.subworkflow(
        name="make_plots_and_report",
        func="single_cell.workflows.cohort_qc.create_cohort_oncoplot",
        args=(
            config,
            mgd.InputFile(cbiofile_paths["filtered_germline_maf"]),
            mgd.InputFile(cbiofile_paths["annotated_somatic_maf"]),
            mgd.InputFile(maftools_filepaths["maftools_cna"]),
            mgd.OutputFile(maftools_filepaths["maftools_maf"]),
            mgd.OutputFile(maftools_filepaths["cohort_oncoplot"])
        ),
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='single_cell.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args['out_dir'],
            list(cbiofile_paths.values()) + list(maftools_filepaths.values()),
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': inpututils.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'cohort_qc'}
        }
    )
    pyp.run(workflow)
