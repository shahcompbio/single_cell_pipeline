import os

import pypeliner
import pypeliner.managed as mgd


def create_patient_workflow(
        pseudobulk_group, mafs, sample_all_snv_csvs,
        mutationreport, merged_maf, high_impact_maf, merged_snvs,
        merged_high_impact_snvs
):
    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1, }
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.transform(
        name='merge_mafs',
        func='single_cell.workflows.qc.tasks.merge_mafs',
        args=(
            mafs,
            mgd.OutputFile(merged_maf),
        ),
        kwargs={"id_colname": True}
    )
    workflow.transform(
        name='filter_merged_maf',
        func='single_cell.workflows.qc.tasks.filter_maf_for_high_impact',
        args=(
            mgd.InputFile(merged_maf),
            mgd.OutputFile(high_impact_maf),
        ),
    )
    workflow.transform(
        name='merge_snvs',
        func='single_cell.workflows.qc.tasks.merge_snvs',
        args=(
            sample_all_snv_csvs,
            mgd.OutputFile(merged_snvs),
        ),
        kwargs={"id_colname": True}
    )
    workflow.transform(
        name='filter_snvs',
        func='single_cell.workflows.qc.tasks.filter_snvs_for_high_impact',
        args=(
            mgd.InputFile(merged_snvs),
            mgd.OutputFile(merged_high_impact_snvs),
        )
    )

    workflow.transform(
        name='mutationreport',
        func='single_cell.workflows.qc.tasks.create_mutation_report',
        args=(
            pseudobulk_group,
            mgd.InputFile(merged_maf),
            mgd.InputFile(high_impact_maf),
            mgd.InputFile(merged_high_impact_snvs),
            mgd.OutputFile(mutationreport),
        ),
    )

    return workflow


def create_sample_level_plots(
        patient, cell_id, library_id, mappability_file,
        strelka_file, museq_file, cosmic_status_file, snpeff_file,
        dbsnp_status_file, trinuc_file, counts_file,
        destruct_breakpoint_annotation, destruct_breakpoint_counts,
        lumpy_breakpoint_annotation, lumpy_breakpoint_evidence,
        haplotype_allele_data, annotation_metrics, hmmcopy_reads,
        hmmcopy_segs, hmmcopy_metrics, alignment_metrics, gc_metrics,
        indel_file, reporthtml, maf, snvs_all_csv, out_dir, config
):

    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1, }

    # vep_reference = config['vep']
    scp_qc_docker = config["docker"]

    prefix = os.path.join(out_dir, patient, cell_id, library_id)

    plots_tar = os.path.join(prefix, "qc_plots.tar")

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.transform(
        name='vcf2maf',
        func='single_cell.workflows.qc.tasks.vcf2maf',
        args=(
            mgd.InputFile(indel_file),
            mgd.OutputFile(maf),
            mgd.TempSpace('vcf2maf_temp'),
            "vep_reference",
        ),
        kwargs=(
            {"docker_image": scp_qc_docker["vcf2maf"]}
        )
    )
    print(out_dir, os.path.isdir(out_dir), maf, os.path.exists(maf))
    workflow.transform(
        name='qc_plots',
        func="single_cell.workflows.qc.scripts.single_cell_qc_plots.qc_plots",
        args=(
            cell_id,
            mgd.InputFile(mappability_file),
            mgd.InputFile(strelka_file),
            mgd.InputFile(museq_file),
            mgd.InputFile(cosmic_status_file),
            mgd.InputFile(snpeff_file),
            mgd.InputFile(dbsnp_status_file),
            mgd.InputFile(trinuc_file),
            mgd.InputFile(counts_file),
            mgd.InputFile(destruct_breakpoint_annotation),
            mgd.InputFile(destruct_breakpoint_counts),
            mgd.InputFile(lumpy_breakpoint_annotation),
            mgd.InputFile(lumpy_breakpoint_evidence),
            mgd.InputFile(haplotype_allele_data),
            mgd.InputFile(annotation_metrics),
            mgd.InputFile(hmmcopy_reads),
            mgd.InputFile(hmmcopy_segs),
            mgd.InputFile(hmmcopy_metrics),
            mgd.InputFile(alignment_metrics),
            mgd.InputFile(gc_metrics),
            library_id,
            prefix,
            mgd.OutputFile(snvs_all_csv),
            mgd.TempSpace("qc_plots"),
            mgd.OutputFile(plots_tar)
        ),
    )

    workflow.transform(
        name='create_main_report',
        func="single_cell.workflows.qc.tasks.sample_level_report",
        args=(

            mgd.InputFile(snvs_all_csv),
            mgd.InputFile(plots_tar),
            mgd.InputFile(maf),
            mgd.OutputFile(reporthtml),
            cell_id + "_" + library_id,
        ),
    )
    print(out_dir, os.path.isdir(out_dir), maf, os.path.exists(maf))
    return workflow
