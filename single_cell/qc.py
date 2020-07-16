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
    jira_ids =  {samplegroup: data["jira_id"] for samplegroup, data in sampledata.items()}
    snv_jira_ids =  {samplegroup: data["snv_genotyping_jira_id"] for samplegroup, data in sampledata.items()}
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
        name='create_sample_level_qc_workflow',
        ctx={'mem':450, 'docker_image': config['docker']['single_cell_pipeline']},
        func="single_cell.workflows.qc.create_sample_level_qc_workflow",
        axes=('sample_id', 'library_id', ),
        args=(
            mgd.InputFile('sample_ids', 'sample_id', 'library_id', fnames=sample_ids),
            mgd.InputFile('library_ids', 'sample_id', 'library_id', fnames=library_ids),
            mgd.InputFile('jira_id', 'sample_id', 'library_id', fnames=jira_ids),
            mgd.InputFile('snv_jira_id', 'sample_id', 'library_id', fnames=snv_jira_ids),
            mgd.OutputFile('sample_level_report_htmls', 'sample_id', 'library_id', fnames=sample_level_report_htmls),
            mgd.OutputFile('mafs', 'sample_id', 'library_id', fnames=mafs),
            mgd.OutputFile('snvs_all', 'sample_id', 'library_id', fnames=snvs_all),
            tmp_dir,
            out_dir,
            mgd.InputFile('outpaths', 'sample_id', 'library_id', fnames=outpaths),

        ),
    )

    # sample_all_snv_csvs = tasks.get_snv_all_csvs(os.path.join(out_dir, pseudobulk_group))
    
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