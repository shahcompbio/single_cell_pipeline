import os
import re

import pypeliner.managed as mgd
from single_cell.utils import inpututils
from single_cell.workflows import align

import pypeliner
import sys

from single_cell.utils import inpututils
from single_cell.workflows.qc import tasks

def qc_workflow(args):
    data = inpututils.load_qc_input(args["input_yaml"])
    config = inpututils.load_config(args)
    config = config["sv_genotyping"]

    out_dir = args["out_dir"]
    tmp_dir = args["tmpdir"]

    pseudobulk_groups = data.keys()
    #inputs
    sampledata = {pseudobulk_group: data[pseudobulk_group] for pseudobulk_group in pseudobulk_groups}

    #outputs
    indelreports = {pseudobulk_group: os.path.join(out_dir, pseudobulk_group, "indelreport.html") for pseudobulk_group in pseudobulk_groups}
    mutationreports = {pseudobulk_group: os.path.join(out_dir, pseudobulk_group, "mutationreport.html") for pseudobulk_group in pseudobulk_groups}
    grouplevelmafs = {pseudobulk_group: os.path.join(out_dir, pseudobulk_group, "grouplevelmaf.maf") for pseudobulk_group in pseudobulk_groups}
    grouplevel_high_impact_mafs = {pseudobulk_group: os.path.join(out_dir, pseudobulk_group, "grouplevel_high_impact_maf.maf") for pseudobulk_group in pseudobulk_groups}
    grouplevel_high_impact_merged_snvs = {pseudobulk_group: os.path.join(out_dir, pseudobulk_group, "grouplevel_high_impact_merged_snvs.csv") for pseudobulk_group in pseudobulk_groups}
    grouplevel_snvs = {pseudobulk_group: os.path.join(out_dir, pseudobulk_group, "grouplevel_snvs.csv") for pseudobulk_group in pseudobulk_groups}

    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config['docker']['single_cell_pipeline']}
    )

    workflow.setobj(
        obj=mgd.OutputChunks('pseudobulk_group'),
        value=list(pseudobulk_groups),
    )

    workflow.subworkflow(
        name='pseudobulk_group_level_qc_workflow',
        ctx={'mem':450, 'docker_image': config['docker']['single_cell_pipeline']},
        func="single_cell.qc.pseudobulk_group_level_qc_workflow",
        axes=('pseudobulk_group',),
        args=(
            mgd.InputFile('pseudobulk_group', 'pseudobulk_group', fnames={psg: psg for psg in pseudobulk_groups}),
            data,
            # mgd.OutputFile('indelreport', 'pseudobulk_group', fnames=indelreports),
            mgd.OutputFile('mutationreport', 'pseudobulk_group', fnames=mutationreports),
            mgd.OutputFile('grouplevelmaf', 'pseudobulk_group', fnames=grouplevelmafs),
            mgd.OutputFile('grouplevel_high_impact_maf', 'pseudobulk_group', fnames=grouplevel_high_impact_mafs),
            mgd.OutputFile('grouplevel_snvs', 'pseudobulk_group', fnames=grouplevel_snvs),
            mgd.OutputFile('grouplevel_high_impact_merged_snvs', 'pseudobulk_group', fnames=grouplevel_high_impact_merged_snvs),
            config,
            out_dir,
            tmp_dir,

        ),
    )
 
    return workflow

def pseudobulk_group_level_qc_workflow(pseudobulk_group, data,  mutationreport, grouplevelmaf, 
                                  grouplevel_high_impact_maf, grouplevel_snvs, grouplevel_high_impact_merged_snvs, config, out_dir, tmp_dir):
    sampledata = data[pseudobulk_group]
    samplegroups = sampledata.keys()

    sample_ids = {samplegroup: samplegroup[0] for samplegroup in samplegroups}
    library_ids = {samplegroup: samplegroup[1] for samplegroup in samplegroups}

    mappability_files =  {samplegroup: data["mappability"] for samplegroup, data in sampledata.items()}
    strelka_files =  {samplegroup: data["strelka"] for samplegroup, data in sampledata.items()}   
    museq_files =  {samplegroup: data["museq"] for samplegroup, data in sampledata.items()}
    cosmic_status_files =  {samplegroup: data["cosmic_status"] for samplegroup, data in sampledata.items()}
    snpeff_files =  {samplegroup: data["snpeff"] for samplegroup, data in sampledata.items()}
    dbsnp_status_files =  {samplegroup: data["dbsnp_status"] for samplegroup, data in sampledata.items()}
    trinuc_files =  {samplegroup: data["trinuc"] for samplegroup, data in sampledata.items()}
    counts_files =  {samplegroup: data["counts"] for samplegroup, data in sampledata.items()}
    breakpoint_counts =  {samplegroup: data["breakpoint_counts"] for samplegroup, data in sampledata.items()}
    breakpoint_annotation =  {samplegroup: data["breakpoint_annotation"] for samplegroup, data in sampledata.items()}
    haplotype_allele_data =  {samplegroup: data["haplotype_allele_data"] for samplegroup, data in sampledata.items()}
    annotation_metrics =  {samplegroup: data["annotation_metrics"] for samplegroup, data in sampledata.items()}
    hmmcopy_reads =  {samplegroup: data["hmmcopy_reads"] for samplegroup, data in sampledata.items()}
    hmmcopy_segs =  {samplegroup: data["hmmcopy_segs"] for samplegroup, data in sampledata.items()}
    hmmcopy_metrics =  {samplegroup: data["hmmcopy_metrics"] for samplegroup, data in sampledata.items()}
    alignment_metrics =  {samplegroup: data["alignment_metrics"] for samplegroup, data in sampledata.items()}
    gc_metrics =  {samplegroup: data["gc_metrics"] for samplegroup, data in sampledata.items()}
    indel_files =  {samplegroup: data["indel_file"] for samplegroup, data in sampledata.items()}

    sample_level_report_htmls =  {samplegroup: os.path.join(out_dir, pseudobulk_group, samplegroup[0], samplegroup[1], "mainreport.html") for samplegroup in samplegroups}
    mafs =  {samplegroup: os.path.join(out_dir, pseudobulk_group, samplegroup[0], samplegroup[1], "samplelevelmaf.maf") for samplegroup in samplegroups}
    snvs_all =  {samplegroup: os.path.join(out_dir, pseudobulk_group, samplegroup[0], samplegroup[1], "snvs_all.csv") for samplegroup in samplegroups}
    outpaths = {samplegroup: os.path.join(pseudobulk_group, samplegroup[0], samplegroup[1]) for samplegroup in samplegroups}

    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config['docker']['single_cell_pipeline']}
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'library_id'),
        value=list(samplegroups),
    )


    workflow.subworkflow(
        name='create_sample_level_plots',
        ctx={'mem':450, 'docker_image': config['docker']['single_cell_pipeline']},
        func="single_cell.workflows.qc.create_sample_level_plots",
        axes=('sample_id', 'library_id', ),
        args=(
            mgd.InputFile('sample_ids', 'sample_id', 'library_id', fnames=sample_ids),
            mgd.InputFile('library_ids', 'sample_id', 'library_id', fnames=library_ids),     
            mgd.InputFile('mappability', 'sample_id', 'library_id', fnames=mappability_files),
            mgd.InputFile('strelka', 'sample_id', 'library_id', fnames=strelka_files),
            mgd.InputFile('museq', 'sample_id', 'library_id', fnames=museq_files),
            mgd.InputFile('cosmic_status', 'sample_id', 'library_id', fnames=cosmic_status_files),
            mgd.InputFile('snpeff', 'sample_id', 'library_id', fnames=snpeff_files),
            mgd.InputFile('dbsnp_status', 'sample_id', 'library_id', fnames=dbsnp_status_files),
            mgd.InputFile('trinuc', 'sample_id', 'library_id', fnames=trinuc_files),
            mgd.InputFile('counts', 'sample_id', 'library_id', fnames=counts_files),
            mgd.InputFile('breakpoint_annotation', 'sample_id', 'library_id', fnames=breakpoint_annotation),
            mgd.InputFile('breakpoint_counts', 'sample_id', 'library_id', fnames=breakpoint_counts),
            mgd.InputFile('haplotype_allele_data', 'sample_id', 'library_id', fnames=haplotype_allele_data),
            mgd.InputFile('annotation_metrics', 'sample_id', 'library_id', fnames=annotation_metrics),
            mgd.InputFile('hmmcopy_reads', 'sample_id', 'library_id', fnames=hmmcopy_reads),
            mgd.InputFile('hmmcopy_segs', 'sample_id', 'library_id', fnames=hmmcopy_segs),
            mgd.InputFile('hmmcopy_metrics', 'sample_id', 'library_id', fnames=hmmcopy_metrics),
            mgd.InputFile('alignment_metrics', 'sample_id', 'library_id', fnames=alignment_metrics),
            mgd.InputFile('gc_metrics', 'sample_id', 'library_id', fnames=gc_metrics),       
            mgd.InputFile('indel_files', 'sample_id', 'library_id', fnames=indel_files),       
            mgd.OutputFile('sample_level_report_htmls', 'sample_id', 'library_id', fnames=sample_level_report_htmls),
            mgd.OutputFile('mafs', 'sample_id', 'library_id', fnames=mafs),
            mgd.OutputFile('snvs_all', 'sample_id', 'library_id', fnames=snvs_all),
            tmp_dir,
            out_dir,
            mgd.InputFile('outpaths', 'sample_id', 'library_id', fnames=outpaths),
        ),
    )
    
    workflow.subworkflow(
        name='create_pseudobulk_group_workflow',
        ctx={'mem':450, 'docker_image': config['docker']['single_cell_pipeline']},
        func="single_cell.workflows.qc.create_pseudobulk_group_workflow",
        args=(
            pseudobulk_group,
            mgd.InputFile("mafs", "sample_id", "library_id", fnames=mafs, axes_origin=[]),
            mgd.InputFile("snvs_all", "sample_id", "library_id", fnames=snvs_all, axes_origin=[]),
            mgd.OutputFile(mutationreport),
            mgd.OutputFile(grouplevelmaf),
            mgd.OutputFile(grouplevel_high_impact_maf),
            mgd.OutputFile(grouplevel_snvs),
            mgd.OutputFile(grouplevel_high_impact_merged_snvs),
            
        ),
    )

    return workflow
 


def qc_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = qc_workflow(args)

    pyp.run(workflow)
