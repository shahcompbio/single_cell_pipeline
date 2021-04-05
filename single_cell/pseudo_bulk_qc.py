import os

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import inpututils


def pseudo_bulk_qc_workflow(args):
    data, patient = inpututils.load_qc_input(args["input_yaml"])
    config = inpututils.load_config(args)
    config = config["qc"]
    out_dir = args["out_dir"]

    mutationreports = os.path.join(out_dir, patient, "mutationreport.html")
    grouplevelmafs = os.path.join(out_dir, patient, "supporting_files", "grouplevelmaf.maf")
    grouplevel_high_impact_mafs = os.path.join(out_dir, patient, "supporting_files", "grouplevel_high_impact_maf.maf")
    grouplevel_high_impact_merged_snvs = os.path.join(out_dir, patient, "supporting_files",  "grouplevel_high_impact_merged_snvs.csv")
    grouplevel_snvs = os.path.join(out_dir, patient, "supporting_files",  "grouplevel_snvs.csv")

    mappability_files = {label: paths["mappability"] for label, paths in data.items()}
    strelka_files = {label: paths["strelka"] for label, paths in data.items()}
    museq_files = {label: paths["museq"] for label, paths in data.items()}
    cosmic_status_files = {label: paths["cosmic_status"] for label, paths in data.items()}
    snpeff_files = {label: paths["snpeff"] for label, paths in data.items()}
    dbsnp_status_files = {label: paths["dbsnp_status"] for label, paths in data.items()}
    trinuc_files = {label: paths["trinuc"] for label, paths in data.items()}
    counts_files = {label: paths["counts"] for label, paths in data.items()}
    breakpoint_counts_files = {label: paths["destruct_breakpoint_counts"]
                               for label, paths in data.items()}
    destruct_breakpoint_annotation_files = {label: paths["destruct_breakpoint_annotation"]
                                            for label, paths in data.items()}
    lumpy_breakpoint_annotation_files = {label: paths["lumpy_breakpoint_annotation"]
                                         for label, paths in data.items()}
    lumpy_breakpoint_evidence_files = {label: paths["lumpy_breakpoint_evidence"]
                                       for label, paths in data.items()}
    haplotype_allele_data_files = {label: paths["haplotype_allele_data"]
                                   for label, paths in data.items()}
    annotation_metrics_files = {label: paths["annotation_metrics"]
                                for label, paths in data.items()}
    hmmcopy_reads_files = {label: paths["hmmcopy_reads"] for label, paths in data.items()}
    hmmcopy_segs_files = {label: paths["hmmcopy_segs"] for label, paths in data.items()}
    hmmcopy_metrics_files = {label: paths["hmmcopy_metrics"] for label, paths in data.items()}
    alignment_metrics_files = {label: paths["alignment_metrics"]
                               for label, paths in data.items()}
    gc_metrics_files = {label: paths["gc_metrics"] for label, paths in data.items()}
    indel_files = {label: paths["indel_file"] for label, paths in data.items()}

    label_dir = os.path.join(out_dir, patient, '{sample_id}', '{library_id}')
    
    sample_level_report_htmls = os.path.join(label_dir,  "mainreport.html")
    sample_level_maf = os.path.join(label_dir, "supporting_files", "samplelevelmaf.maf")
    snvs_all = os.path.join(label_dir, "supporting_files", 'snvs_all.csv')

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'library_id', ),
        value=list(data.keys()),
    )

    workflow.subworkflow(
        name='create_sample_level_plots',
        func="single_cell.workflows.pseudo_bulk_qc.create_sample_level_plots",
        axes=('sample_id', 'library_id',),
        args=(
            patient,
            mgd.InputInstance('sample_id'),
            mgd.InputInstance('library_id'),
            mgd.InputFile('mappability', 'sample_id', 'library_id',
                          fnames=mappability_files),
            mgd.InputFile('strelka', 'sample_id', 'library_id', fnames=strelka_files),
            mgd.InputFile('museq', 'sample_id', 'library_id', fnames=museq_files),
            mgd.InputFile('cosmic_status', 'sample_id', 'library_id',
                          fnames=cosmic_status_files),
            mgd.InputFile('snpeff', 'sample_id', 'library_id', fnames=snpeff_files),
            mgd.InputFile('dbsnp_status', 'sample_id', 'library_id',
                          fnames=dbsnp_status_files),
            mgd.InputFile('trinuc', 'sample_id', 'library_id', fnames=trinuc_files),
            mgd.InputFile('counts', 'sample_id', 'library_id', fnames=counts_files),
            mgd.InputFile('destruct_breakpoint_annotation', 'sample_id', 'library_id',
                          fnames=destruct_breakpoint_annotation_files),
            mgd.InputFile('destruct_breakpoint_counts', 'sample_id', 'library_id',
                          fnames=breakpoint_counts_files),
            mgd.InputFile('lumpy_breakpoint_annotation', 'sample_id', 'library_id',
                          fnames=lumpy_breakpoint_annotation_files),
            mgd.InputFile('lumpy_breakpoint_evidence', 'sample_id', 'library_id',
                          fnames=lumpy_breakpoint_evidence_files),
            mgd.InputFile('haplotype_allele_data', 'sample_id', 'library_id',
                          fnames=haplotype_allele_data_files),
            mgd.InputFile('annotation_metrics', 'sample_id', 'library_id',
                          fnames=annotation_metrics_files),
            mgd.InputFile('hmmcopy_reads', 'sample_id', 'library_id', fnames=hmmcopy_reads_files),
            mgd.InputInstance("sample_id"),
            mgd.InputFile('hmmcopy_segs', 'sample_id', 'library_id', fnames=hmmcopy_segs_files),
            mgd.InputFile('hmmcopy_metrics', 'sample_id', 'library_id', fnames=hmmcopy_metrics_files),
            mgd.InputFile('alignment_metrics', 'sample_id', 'library_id',
                          fnames=alignment_metrics_files),
            mgd.InputFile('gc_metrics', 'sample_id', 'library_id', fnames=gc_metrics_files),
            mgd.InputFile('indel_files', 'sample_id', 'library_id', fnames=indel_files),
            mgd.OutputFile('sample_level_report_htmls', 'sample_id', 'library_id',
                           template=sample_level_report_htmls),
            mgd.OutputFile('mafs', 'sample_id', 'library_id', template=sample_level_maf),
            mgd.OutputFile('snvs_all', 'sample_id', 'library_id', template=snvs_all),
            out_dir,
            config

        ),
    )
    workflow.subworkflow(
        name='create_patient_workflow',
        func="single_cell.workflows.pseudo_bulk_qc.create_patient_workflow",
        args=(
            patient,
            mgd.InputFile("mafs", "sample_id", "library_id",
                          template=sample_level_maf, axes_origin=[]),
            mgd.InputFile("snvs_all", "sample_id", "library_id",
                          template=snvs_all, axes_origin=[]),
            mutationreports,
            grouplevelmafs,
            grouplevel_high_impact_mafs,
            grouplevel_snvs,
            grouplevel_high_impact_merged_snvs,
        ),
    )

    return workflow


def pseudo_bulk_qc_pipeline(args):
    args["out_dir"] = os.path.realpath(args["out_dir"])
    args["tmpdir"] = os.path.realpath(args["tmpdir"])

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = pseudo_bulk_qc_workflow(args)

    pyp.run(workflow)
