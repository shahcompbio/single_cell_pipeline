import pypeliner
import pypeliner.managed as mgd


def cna_annotation_workflow(
        hmmcopy_dict,
        output_cbio_table,
        output_maftools_table,
        output_segs,
        gtf
):
    """Run classify copynumber on hmmcopy dictionary and generate a cohort cna table and file containg segments.

    Args:
        config ([dict]): [config]
        hmmcopy_dict ([dict]): [dictionary of sample: hmmcopy file]
        output_cbio_table ([str]): [path to cna data for cbio]
        output_maftools_table ([str]): [path to cna data for maftools]
        output_segs ([str]): [path to output segments for cbio]
        gtf ([str]): [path to gtf file]

    Returns:
    """
    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_label', 'library_label'),
        value=list(hmmcopy_dict.keys()),
    )

    workflow.transform(
        name='classify_hmmcopy_files',
        func='single_cell.workflows.cohort_qc.tasks.classify_hmmcopy',
        axes=("sample_label",),
        args=(
            mgd.InputInstance("sample_label"),
            mgd.InputFile(
                'hmmcopy', 'sample_label', 'library_label', fnames=hmmcopy_dict, axes_origin=[]
            ),
            gtf,
            mgd.TempSpace("annotated_maf_tmp", "sample_label"),
            mgd.TempOutputFile('amps', 'sample_label'),
            mgd.TempOutputFile('dels', 'sample_label'),
        ),
    )

    workflow.transform(
        name='merge_amp_tables',
        func='single_cell.workflows.cohort_qc.tasks.merge_cna_tables',
        args=(
            mgd.TempInputFile('amps', 'sample_label', axes_origin=[]),
            mgd.TempOutputFile("merged_amps"),
        ),
    )

    workflow.transform(
        name='merge_del_tables',
        func='single_cell.workflows.cohort_qc.tasks.merge_cna_tables',
        args=(
            mgd.TempInputFile('dels', 'sample_label', axes_origin=[]),
            mgd.TempOutputFile("merged_dels"),
        ),
    )

    workflow.transform(
        name='make_cbio_cna_table',
        func='single_cell.workflows.cohort_qc.tasks.make_cbio_cna_table',
        args=(
            mgd.TempInputFile('merged_amps'),
            mgd.TempInputFile('merged_dels'),
            mgd.OutputFile(output_cbio_table),
        ),
    )

    workflow.transform(
        name='make_maftools_cna_table',
        func='single_cell.workflows.cohort_qc.tasks.make_maftools_cna_table',
        args=(
            mgd.TempInputFile('merged_amps'),
            mgd.TempInputFile('merged_dels'),
            output_maftools_table
        ),
    )

    workflow.transform(
        name='generate_segmental_copynumber',
        func='single_cell.workflows.cohort_qc.tasks.generate_segmental_copynumber',
        axes=("sample_label",),
        args=(
            mgd.InputFile('hmmcopy', 'sample_label', 'library_label', fnames=hmmcopy_dict),
            mgd.TempInputObj(
                'sample_ids', 'sample_label', 'library_label',
            ),
            mgd.TempOutputFile('segmental_cn', 'sample_label'),
            mgd.InputInstance('sample_label'),
        ),
    )

    workflow.transform(
        name='merge_segmental_cn',
        func='single_cell.workflows.cohort_qc.tasks.merge_segmental_cn',
        args=(
            mgd.TempInputFile('segmental_cn', 'sample_label', axes_origin=[]),
            mgd.OutputFile(output_segs)
        ),
    )

    return workflow


def preprocess_mafs_workflow(
        germline_mafs, somatic_mafs,
        cohort_germline_maf, cohort_somatic_maf, api_key
):
    """Take germline and somatic mafs, annotates, filters and merges them.

    Args:
        config ([dict]): [config]
        germline_mafs ([dic]): [sample: germline maf path]
        somatic_mafs ([dict]): [sample: somatic maf path]
        cohort_germline_maf ([str]): [merged germline output]
        cohort_somatic_maf ([str]): [merged somatic output]
        api_key ([str]): [API key for onckb annotator]

    Returns:
    """
    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_label', ),
        value=list(somatic_mafs.keys()),
    )

    workflow.transform(
        name='annotate_germline_mafs',
        func='single_cell.workflows.cohort_qc.tasks.annotate_maf_with_oncokb',
        axes=("sample_label",),
        args=(
            mgd.InputFile(
                'germlne_maf', 'sample_label', fnames=germline_mafs
            ),
            api_key,
            mgd.TempSpace("annotated_germline_maf_tmp", 'sample_label'),
            mgd.TempOutputFile("annotated_germline_maf", 'sample_label'),
        ),
    )

    workflow.transform(
        name='filter_germline_mafs',
        func='single_cell.workflows.cohort_qc.tasks.filter_maf',
        axes=("sample_label",),
        args=(
            mgd.TempInputFile("annotated_germline_maf", 'sample_label'),
            mgd.TempOutputFile("filtered_germline_maf", "sample_label")
        ),
    )

    workflow.transform(
        name='merge_germline_mafs',
        func='single_cell.workflows.cohort_qc.tasks.merge_mafs',
        args=(
            mgd.TempInputFile(
                "filtered_germline_maf", 'sample_label', axes_origin=[]
            ),
            mgd.OutputFile(cohort_germline_maf)
        ),
    )

    workflow.transform(
        name='annotate_somatic_mafs',
        func='single_cell.workflows.cohort_qc.tasks.annotate_maf_with_oncokb',
        axes=("sample_label",),
        args=(
            mgd.InputFile('somatic_maf', 'sample_label', fnames=somatic_mafs),
            api_key,
            mgd.TempSpace("annotated_somatic_maf_tmp", 'sample_label'),
            mgd.TempOutputFile("annotated_somatic_maf", 'sample_label'),
        ),
    )
    workflow.transform(
        name='merge_somatic_mafs',
        func='single_cell.workflows.cohort_qc.tasks.merge_mafs',
        args=(
            mgd.TempInputFile(
                "annotated_somatic_maf", 'sample_label', axes_origin=[]
            ),
            mgd.OutputFile(cohort_somatic_maf)
        ),
    )

    return workflow


def create_cohort_oncoplot(
        config, merged_germline,
        merged_somatic, maftools_cna, maftools_maf, oncoplot
):
    """create oncoplot from cna table, and germlinne/somatic dataa.

    Args:
        config ([dict]): [config]
        merged_germline ([str]): [path to merged germline file]
        merged_somatic ([str]): [path to merged somatic file]
        maftools_cna ([str]): [path to merged cna data]
        maftools_maf ([sstr]): [path to output prepped maftools input maf]
        oncoplot ([str]): [path to output oncoplot]

    Returns:
        [type]: [description]
    """
    workflow = pypeliner.workflow.Workflow()

    non_synonymous_labels = config["non_synonymous_labels"]

    workflow.transform(
        name='postprocess_maf',
        func='single_cell.workflows.cohort_qc.tasks.prepare_maf_for_maftools',
        args=(
            mgd.InputFile(merged_germline),
            mgd.InputFile(merged_somatic),
            mgd.OutputFile(maftools_maf),
            non_synonymous_labels,
            mgd.TempOutputFile("vcNames")
        ),
    )

    workflow.transform(
        name='make_oncoplot',
        func='single_cell.workflows.cohort_qc.tasks.make_oncoplot',
        args=(
            mgd.InputFile(maftools_maf),
            mgd.InputFile(maftools_cna),
            mgd.OutputFile(oncoplot),
            mgd.TempInputFile("vcNames")

        ),
    )

    return workflow
