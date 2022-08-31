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
    report = os.path.join(
        root_dir, "maftools",  "report.html"
    )
    return {
        "cohort_oncoplot": cohort_oncoplot,
        "maftools_maf": maftools_maf,
        "maftools_cna": maftools_cna,
        "report": report
    }


def cohort_qc_pipeline(args):
    """Process maf, run classify copynumber, make plots.
    Args:
        args ([dict]): [pipeline arguments]
    """
    config = inpututils.load_config(args)
    config = config["cohort_qc"]

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = pypeliner.workflow.Workflow()

    out_dir = args["output_prefix"]
    api_key = args["API_key"]

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    # inputs
    cohort, germline_mafs, vcfs, hmmcopy = inpututils.load_cohort_qc_inputs(
        args["input_yaml"]
    )

    museq = {
        label: data["museq"] for label, data in vcfs.items()
    }
    strelka_snv = {
        label: data["strelka_snv"] for label, data in vcfs.items()
    }
    strelka_indel = {
        label: data["strelka_indel"] for label, data in vcfs.items()
    }
    hmmcopy_files = {
        label: data["hmmcopy"] for label, data in hmmcopy.items()
    }
    hmmcopy_metrics_files = {
        label: data["hmmcopy_metrics"] for label, data in hmmcopy.items()
    }
    # outputs
    cbiofile_paths = get_cbioportal_paths(os.path.join(out_dir, cohort))
    maftools_filepaths = get_maftools_paths(os.path.join(out_dir, cohort))

    workflow.setobj(
        obj=mgd.OutputChunks('sample_label', 'library_label'),
        value=list(museq.keys()),
    )
    workflow.subworkflow(
        name="merge_somatic_mafs",
        func="single_cell.workflows.cohort_qc.merge_somatic_mafs",
        axes=('sample_label',),
        args=(
            mgd.InputInstance('sample_label'),
            config,
            mgd.InputFile(
                'museq', 'sample_label', 'library_label',
                fnames=museq, axes_origin=[]
            ),
            mgd.InputFile(
                'strelka_snv', 'sample_label', 'library_label',
                fnames=strelka_snv, axes_origin=[]
            ),
            mgd.InputFile(
                'strelka_indel', 'sample_label', 'library_label',
                fnames=strelka_indel, axes_origin=[]
            ),
            mgd.TempOutputFile('somatic_maf', 'sample_label')
        ),
    )
    
    workflow.subworkflow(
        name="classifycopynumber",
        func="single_cell.workflows.cohort_qc.cna_annotation_workflow",
        args=(
            config,
            mgd.InputFile(
                'hmmcopy_dict', 'sample_label', 'library_label',
                fnames=hmmcopy_files, axes_origin=[]
            ),
            mgd.InputFile(
                'hmmcopy_metrics_dict', 'sample_label', 'library_label',
                fnames=hmmcopy_metrics_files, axes_origin=[]
            ),
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
            mgd.TempInputFile(
                'somatic_maf',  'sample_label',
                axes_origin=[]
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
            mgd.OutputFile(maftools_filepaths["cohort_oncoplot"]),
            mgd.OutputFile(maftools_filepaths["report"]),
            cohort
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
