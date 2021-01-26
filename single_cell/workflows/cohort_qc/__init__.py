import logging
import os

import pypeliner
import pypeliner.managed as mgd
from wgs.config import config


def cna_annotation_workflow(config, hmmcopy_dict, output_cbio_table, 
    output_maftools_table, output_segs, gtf):
    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config['docker']['single_cell_pipeline']}
    )

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
            mgd.InputFile('hmmcopy', 'sample_label', 'library_label', fnames=hmmcopy_dict, 
                axes_origin=[]),
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
            mgd.TempOutputFile('segmental_cn', 'sample_label'),
            mgd.InputInstance('sample_label')
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


def preprocess_mafs_workflow(config, germline_mafs, somatic_mafs, cohort_maf, api_key):
    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config['docker']['single_cell_pipeline']}
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_label',),
        value=list(somatic_mafs.keys()),
    )

    workflow.transform(
        name='annotate_germline_mafs',
        func='single_cell.workflows.cohort_qc.tasks.annotate_maf_with_oncokb',
        axes=("sample_label",),
        args=(
            mgd.InputFile('germlne_maf', 'sample_label', fnames=germline_mafs),
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
        name='annotate_germline_class',
        func='single_cell.workflows.cohort_qc.tasks.annotate_germline_somatic',
        axes=("sample_label",),
        args=(
            mgd.TempInputFile('filtered_germline_maf', 'sample_label'),
            mgd.TempOutputFile("filtered_class_labeled_germline_maf", 'sample_label'),
            True
        ),
    )

    workflow.transform(
        name='annotate_somatic_class',
        func='single_cell.workflows.cohort_qc.tasks.annotate_germline_somatic',
        axes=("sample_label",),
        args=(
            mgd.TempInputFile('annotated_somatic_maf', 'sample_label'),
            mgd.TempOutputFile("class_labeled_somatic_maf", 'sample_label'),
            False
        ),
    )

    workflow.transform(
        name='merge_filtered_germline_somatic',
        func='single_cell.workflows.cohort_qc.tasks.merge_mafs',
        args=(
            mgd.TempInputFile("filtered_class_labeled_germline_maf", 'sample_label', axes_origin=[]),
            mgd.TempInputFile('class_labeled_somatic_maf', 'sample_label', axes_origin=[]),
            mgd.OutputFile(cohort_maf),
        ),
    )

    return workflow


def create_cohort_oncoplot(config, cohort, out_dir,  cohort_maf, cna_table, oncoplot):

    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config['docker']['single_cell_pipeline']}
    )

    non_synonymous_labels = config["non_synonymous_labels"]

    workflow.transform(
        name='postprocess_maf',
        func='single_cell.workflows.cohort_qc.tasks.prepare_maf_for_maftools',
        args=(
            cohort,
            mgd.InputFile(cohort_maf),
            mgd.TempOutputFile("prepared_maf"),
            non_synonymous_labels,
            mgd.TempOutputFile("vcNames")
        ),
    )
    
    workflow.transform(
        name='make_oncoplot',
        func='single_cell.workflows.cohort_qc.tasks.make_oncoplot',
        args=(
            mgd.TempInputFile("prepared_maf"),
            mgd.InputFile(cna_table),
            mgd.OutputFile(oncoplot),
            mgd.TempInputFile("vcNames")

        ),
        kwargs={'docker_image':config["docker"]["pseudo_bulk_qc_html_report"] },
    )

    return workflow