import logging
import os

import pypeliner
import pypeliner.managed as mgd
from wgs.config import config


# def classify_hmmcopy(sample, hmmcopy_dict, gtf, cna_table):
#     workflow = pypeliner.workflow.Workflow(
#         ctx={'docker_image': config.containers('wgs')}
#     )
#     workflow.setobj(
#         obj=mgd.OutputChunks('library_label'),
#         value=list(hmmcopy_dict.keys()),
#     )

#     workflow.transform(
#         name='classify_hmmcopy',
#         func='single_cell.workflows.cohort_qc.tasks.classify_hmmcopy',
#         axes=("library_label",),
#         args=(
#             mgd.InputInstance('library_label'),
#             mgd.InputFile('remixt', 'library_label', fnames=hmmcopy_dict),
#             gtf,
#             mgd.TempSpace('annotated_maf_tmp', 'library_label'),
#             mgd.TempOutputFile('amps', 'library_label'),
#             mgd.TempOutputFile('dels', 'library_label'),
#         ),
#     )

#     workflow.transform(
#         name='make_cna_tables',
#         func='single_cell.workflows.cohort_qc.tasks.make_cna_table',
#         args=(
#             sample,
#             mgd.TempInputFile('amps', 'library_label', axes_origin=[]),
#             mgd.TempInputFile('dels', 'library_label', axes_origin=[]),
#             mgd.OutputFile(cna_table, extensions=['.yaml']),
#         ),
#     )

#     return workflow


def cna_annotation_workflow(config, hmmcopy_dict, output_table, output_segs, gtf):
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
            mgd.InputFile('hmmcopy', 'sample_label', 'library_label', fnames=hmmcopy_dict, axes_origin=[]),
            gtf,
            mgd.TempSpace("annotated_maf_tmp", "sample_label"),
            mgd.TempOutputFile('amps', 'sample_label'),
            mgd.TempOutputFile('dels', 'sample_label'),   
        ),     
    )

    workflow.transform(
        name='make_cna_tables',
        func='single_cell.workflows.cohort_qc.tasks.make_cna_table',
        args=(
            mgd.TempInputFile('amps', 'sample_label', axes_origin=[]),
            mgd.TempInputFile('dels', 'sample_label', axes_origin=[]),
            mgd.OutputFile(output_table, extensins=['.yaml']),
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
        value=list(germline_mafs.keys()),
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
            mgd.TempInputFile('filtered_class_labeled_germline_maf', 'sample_label', axes_origin=[]),
            mgd.TempInputFile('class_labeled_somatic_maf', 'sample_label', axes_origin=[]),
            mgd.OutputFile(cohort_maf),
        ),
    )

    return workflow


def create_cohort_oncoplot(config, cohort, out_dir,  cohort_maf, cna_table, oncoplot):

    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config['docker']['single_cell_pipeline']}
    )


    non_synonymous_labels=["Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", 
        "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation", 
        "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"
    ]

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
        name='format_cna_table',
        func='single_cell.workflows.cohort_qc.tasks.format_cna_table',
        args=(
            mgd.InputFile(cna_table),
            mgd.TempOutputFile("cna_table_formatted")
        ),
    )
    
    workflow.transform(
        name='make_oncoplot',
        func='single_cell.workflows.cohort_qc.tasks.make_oncoplot',
        args=(
            mgd.TempInputFile("prepared_maf"),
            mgd.TempInputFile("cna_table_formatted"),
            mgd.OutputFile(oncoplot),
            mgd.TempInputFile("vcNames")

        ),
        # kwargs={'docker_image':config.containers("wgs_qc_html") },
    )

    return workflow