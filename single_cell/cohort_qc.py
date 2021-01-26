import os
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import inpututils


def cohort_qc_pipeline(args):

    config = inpututils.load_config(args)
    config = config["cohort_qc"]

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config['docker']['single_cell_pipeline']}
    )

    cohort, mafs, hmmcopy = inpututils.load_cohort_qc_inputs(args["input_yaml"])
    out_dir = args["out_dir"]
    tmp_dir = args["tmpdir"]
    api_key = args["API_key"]
    gtf = config["gtf"]

    germline_mafs = {label: data["germline_maf"] for label, data in mafs.items()}
    somatic_mafs = {label: data["somatic_maf"] for label, data in mafs.items()}

    hmmcopy_files = {label: data["hmmcopy"] for label, data in hmmcopy.items()}

    cna_cbioportal_table = os.path.join(out_dir, cohort, "cna_table.tsv.gz")
    segments = os.path.join(out_dir, cohort, "segments.tsv.gz")

    cohort_maf_oncogenic_filtered = os.path.join(out_dir, cohort, "cohort_oncogenic_filtered.maf")
    cohort_oncoplot = os.path.join(out_dir, cohort, "cohort_oncoplot.png")

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
            mgd.InputFile(cna),
            mgd.OutputFile(cohort_oncoplot)
        ),
    )   

    pyp.run(workflow)
