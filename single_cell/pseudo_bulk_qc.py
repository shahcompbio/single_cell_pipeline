import os
import re
import pypeliner.managed as mgd
from single_cell.workflows import align
import pypeliner
import sys
from single_cell.utils import inpututils
from single_cell.workflows.qc import tasks

def qc_workflow(args):
    data = inpututils.load_qc_input(args["input_yaml"])
    config = inpututils.load_config(args)
    config = config["qc"]

    out_dir = args["out_dir"]

    patients = [k[0] for k in data.keys()]

    #outputs
    indelreports = {patient: os.path.join(out_dir, patient, "indelreport.html") 
        for patient in patients
    }
    mutationreports = {patient: os.path.join(out_dir, patient, "mutationreport.html") 
        for patient in patients
    }
    grouplevelmafs = {patient: os.path.join(out_dir, patient, "grouplevelmaf.maf") 
        for patient in patients
    }
    grouplevel_high_impact_mafs = {patient: os.path.join(out_dir, patient, 
        "grouplevel_high_impact_maf.maf") 
        for patient in patients
    }
    grouplevel_high_impact_merged_snvs = {patient: os.path.join(out_dir, patient, 
        "grouplevel_high_impact_merged_snvs.csv") for patient in patients
    }
    grouplevel_snvs = {patient: os.path.join(out_dir, patient, "grouplevel_snvs.csv") 
        for patient in patients
    }

    mappability_files = {label:paths["mappability"] for label, paths in data.items()}
    strelka_files = {label:paths["strelka"] for label, paths in data.items()}
    museq_files = {label:paths["museq"] for label, paths in data.items()}
    cosmic_status_files = {label:paths["cosmic_status"] for label, paths in data.items()}
    snpeff_files = {label:paths["snpeff"] for label, paths in data.items()}
    dbsnp_status_files = {label:paths["dbsnp_status"] for label, paths in data.items()}
    trinuc_files = {label:paths["trinuc"] for label, paths in data.items()}
    counts_files = {label:paths["counts"] for label, paths in data.items()}
    breakpoint_counts_files = {label:paths["destruct_breakpoint_counts"] 
        for label, paths in data.items()}
    destruct_breakpoint_annotation_files = {label:paths["destruct_breakpoint_annotation"] 
        for label, paths in data.items()}
    lumpy_breakpoint_annotation_files = {label:paths["lumpy_breakpoint_annotation"] 
        for label, paths in data.items()}
    lumpy_breakpoint_evidence_files = {label:paths["lumpy_breakpoint_evidence"] 
        for label, paths in data.items()}
    haplotype_allele_data_files = {label:paths["haplotype_allele_data"] 
        for label, paths in data.items()}
    annotation_metrics_files = {label:paths["annotation_metrics"] 
        for label, paths in data.items()}
    hmmcopy_reads_files = {label:paths["hmmcopy_reads"] for label, paths in data.items()}
    hmmcopy_segs_files = {label:paths["hmmcopy_segs"] for label, paths in data.items()}
    hmmcopy_metrics_files = {label:paths["hmmcopy_metrics"] for label, paths in data.items()}
    alignment_metrics_files = {label:paths["alignment_metrics"] 
        for label, paths in data.items()}
    gc_metrics_files = {label:paths["gc_metrics"] for label, paths in data.items()}
    indel_files = {label:paths["indel_file"] for label, paths in data.items()}

    sample_level_report_htmls = {label:os.path.join(out_dir, *label,  "mainreport.html") 
        for label, paths in data.items()}
    sample_level_maf = {label:os.path.join(out_dir, *label,  "samplelevelmaf.maf") 
        for label, paths in data.items()}
    snvs_all = {label:os.path.join(out_dir, *label,  "snvs_all.csv") for label, paths in data.items()}


    patient_names = {patient: patient for patient in patients}

    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config['docker']['single_cell_pipeline']}
    )

    workflow.setobj(
        obj=mgd.OutputChunks('patient', 'sample_id', 'library_id',),
        value=list(data.keys()),
    )    

    workflow.subworkflow(
        name='create_sample_level_plots',
        func="single_cell.workflows.qc.create_sample_level_plots",
        axes=('patient', 'sample_id', 'library_id',),
        args=(
            mgd.InputInstance('patient'),
            mgd.InputInstance('sample_id'),
            mgd.InputInstance('library_id'),
            mgd.InputFile('mappability', 'patient', 'sample_id', 'library_id', 
                fnames=mappability_files),
            mgd.InputFile('strelka', 'patient', 'sample_id', 'library_id', fnames=strelka_files),
            mgd.InputFile('museq', 'patient', 'sample_id', 'library_id', fnames=museq_files),
            mgd.InputFile('cosmic_status', 'patient', 'sample_id', 'library_id', 
                fnames=cosmic_status_files),
            mgd.InputFile('snpeff', 'patient', 'sample_id', 'library_id', fnames=snpeff_files),
            mgd.InputFile('dbsnp_status', 'patient', 'sample_id', 'library_id', 
                fnames=dbsnp_status_files),
            mgd.InputFile('trinuc', 'patient', 'sample_id', 'library_id', fnames=trinuc_files),
            mgd.InputFile('counts', 'patient', 'sample_id', 'library_id', fnames=counts_files),
            mgd.InputFile('destruct_breakpoint_annotation', 'patient', 'sample_id', 'library_id', 
                fnames=destruct_breakpoint_annotation_files),
            mgd.InputFile('destruct_breakpoint_counts','patient',  'sample_id', 'library_id', 
                fnames=breakpoint_counts_files),
            mgd.InputFile('lumpy_breakpoint_annotation', 'patient', 'sample_id', 'library_id', 
                fnames=lumpy_breakpoint_annotation_files),
            mgd.InputFile('lumpy_breakpoint_evidence','patient',  'sample_id', 'library_id', 
                fnames=lumpy_breakpoint_evidence_files),
            mgd.InputFile('haplotype_allele_data', 'patient', 'sample_id', 'library_id', 
                fnames=haplotype_allele_data_files),
            mgd.InputFile('annotation_metrics', 'patient', 'sample_id', 'library_id', 
                fnames=annotation_metrics_files),
            mgd.InputFile('hmmcopy_reads', 'patient', 'sample_id', 'library_id', fnames=hmmcopy_reads_files),
            mgd.InputFile('hmmcopy_segs', 'patient', 'sample_id', 'library_id', fnames=hmmcopy_segs_files),
            mgd.InputFile('hmmcopy_metrics', 'patient', 'sample_id', 'library_id', fnames=hmmcopy_metrics_files),
            mgd.InputFile('alignment_metrics', 'patient', 'sample_id', 'library_id', 
                fnames=alignment_metrics_files),
            mgd.InputFile('gc_metrics','patient',  'sample_id', 'library_id', fnames=gc_metrics_files),       
            mgd.InputFile('indel_files','patient', 'sample_id', 'library_id', fnames=indel_files),  
            mgd.OutputFile('sample_level_report_htmls', 'patient', 'sample_id', 'library_id', 
                fnames=sample_level_report_htmls),
            mgd.OutputFile('mafs', 'patient', 'sample_id', 'library_id', fnames=sample_level_maf),
            mgd.OutputFile('snvs_all', 'patient', 'sample_id', 'library_id',fnames=snvs_all),
            out_dir,
            config  


        ),
    )
    workflow.subworkflow(
        name='create_patient_workflow',
        func="single_cell.workflows.qc.create_patient_workflow",
        axes=('patient',),
        args=(
            mgd.InputInstance('patient'),
            mgd.InputFile("mafs", "patient", "sample_id", "library_id", 
                fnames=sample_level_maf, axes_origin=[]),
            mgd.InputFile("snvs_all", "patient", "sample_id", "library_id", 
                fnames=snvs_all, axes_origin=[]),
            mgd.OutputFile('mutationreport', 'patient', fnames=mutationreports),
            mgd.OutputFile('grouplevelmaf', 'patient', fnames=grouplevelmafs),
            mgd.OutputFile('grouplevel_high_impact_maf', 'patient', 
                fnames=grouplevel_high_impact_mafs
            ),
            mgd.OutputFile('grouplevel_snvs', 'patient', fnames=grouplevel_snvs),
            mgd.OutputFile('grouplevel_high_impact_merged_snvs', 'patient', 
                fnames=grouplevel_high_impact_merged_snvs
            ),
        ),
    )

    return workflow
    
    

def qc_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = qc_workflow(args)

    pyp.run(workflow)


