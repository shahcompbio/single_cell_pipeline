# Douglas Abrams
# 10/14/19

import pypeliner
import pypeliner.managed as mgd
import logging

def create_sv_genotyping_workflow(
        reference,
        input_bam,
        lumpy_csv,
        destruct_csv,
        ao_output,
        ap_output,
        as_output,
        asc_output,
        dp_output,
        gq_output,
        qa_output,
        qr_output,
        ro_output,
        rp_output,
        rs_output,
        sq_output,
        gl_output,
        ab_output,
        config
):
    img = config["docker"]["single_cell_pipeline"]
    genotyping_ctx = {'docker_image': img}

    logging.getLogger("SV_GENOTYPING").error("EXPERIMENTAL FEATURE")

    workflow = pypeliner.workflow.Workflow(ctx=genotyping_ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=list(input_bam.keys()),
    )

    workflow.transform(
        name='lumpy_csv_to_vcf',
        func='single_cell.workflows.sv_genotyping.tasks.varcalls_to_svtyper_input',
        args=(
            mgd.InputFile(lumpy_csv),
            mgd.TempOutputFile("lumpy.vcf"),
            mgd.TempSpace("lumpy_input"),
            "lumpy",
        )
    )

    workflow.transform(
        name='genotype-lumpy',
        func='single_cell.workflows.sv_genotyping.tasks.genotype',
        axes=('cell_id',),
        args=(
            mgd.InputFile("cell_specific.bam", "cell_id", fnames=input_bam, extensions=['.bai']),
            reference,
            mgd.TempInputFile("lumpy.vcf"),
            mgd.TempOutputFile("genotypedlumpy.vcf", "cell_id"),
            mgd.TempOutputFile("genotypedlumpy.csv.gz", "cell_id", extensions=['.yaml']),
            mgd.TempSpace("temp", 'cell_id'),
            mgd.InputInstance('cell_id'),
            config["docker"]["svtyper"]
        )
    )

    workflow.transform(
        name='merge_cell_specific_csv-lumpy',
        func="single_cell.workflows.sv_genotyping.tasks.merge_csvs",
        args=(
            mgd.TempInputFile("genotypedlumpy.csv.gz", 'cell_id', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile("genotypedmerged_lumpy.csv.gz", extensions=['.yaml']),
        )
    )

    workflow.transform(
        name='destruct_csv_to_vcf',
        func='single_cell.workflows.sv_genotyping.tasks.varcalls_to_svtyper_input',
        args=(
            mgd.InputFile(destruct_csv),
            mgd.TempOutputFile("destruct.vcf"),
            mgd.TempSpace("destruct_input"),
            "destruct",
        )
    )

    workflow.transform(
        name='genotype-destruct',
        func='single_cell.workflows.sv_genotyping.tasks.genotype',
        axes=('cell_id',),
        args=(
            mgd.InputFile("cell_specific.bam", "cell_id", fnames=input_bam, extensions=['.bai']),
            reference,
            mgd.TempInputFile("destruct.vcf"),
            mgd.TempOutputFile("genotypeddestruct.vcf", "cell_id"),
            mgd.TempOutputFile("genotypeddestruct.csv.gz", "cell_id", extensions=['.yaml']),
            mgd.TempSpace("tempdestruct", 'cell_id'),
            mgd.InputInstance('cell_id'),
            config["docker"]["svtyper"]
        )
    )

    workflow.transform(
        name='merge_cell_specific_csv-destruct',
        func="single_cell.workflows.sv_genotyping.tasks.merge_csvs",
        args=(
            mgd.TempInputFile("genotypeddestruct.csv.gz", 'cell_id', axes_origin=[], extensions=['.yaml']),
            mgd.TempOutputFile("genotypedmerged_destruct.csv.gz", extensions=['.yaml']),
        )
    )

    workflow.transform(
        name='merge_lumpy_destruct_vcfs',
        func="single_cell.workflows.sv_genotyping.tasks.merge_csvs",
        args=(
            [
                mgd.TempInputFile("genotypedmerged_lumpy.csv.gz", extensions=['.yaml']),
                mgd.TempInputFile("genotypedmerged_destruct.csv.gz", extensions=['.yaml'])
            ],
            mgd.TempOutputFile("lumpy_destruct_merged.csv.gz", extensions=['.yaml']),
        )
    )

    workflow.transform(
        name='write_svtyper_annotations',
        func='single_cell.workflows.sv_genotyping.tasks.write_svtyper_annotations',
        args=(
            mgd.TempInputFile("lumpy_destruct_merged.csv.gz", extensions=['.yaml']),
            {
                'AO': mgd.OutputFile(ao_output, extensions=['.yaml']),
                'AP': mgd.OutputFile(ap_output, extensions=['.yaml']),
                'AS': mgd.OutputFile(as_output, extensions=['.yaml']),
                'ASC': mgd.OutputFile(asc_output, extensions=['.yaml']),
                'DP': mgd.OutputFile(dp_output, extensions=['.yaml']),
                'GQ': mgd.OutputFile(gq_output, extensions=['.yaml']),
                'QA': mgd.OutputFile(qa_output, extensions=['.yaml']),
                'QR': mgd.OutputFile(qr_output, extensions=['.yaml']),
                'RO': mgd.OutputFile(ro_output, extensions=['.yaml']),
                'RP': mgd.OutputFile(rp_output, extensions=['.yaml']),
                'RS': mgd.OutputFile(rs_output, extensions=['.yaml']),
                'SQ': mgd.OutputFile(sq_output, extensions=['.yaml']),
                'GL': mgd.OutputFile(gl_output, extensions=['.yaml']),
                'AB': mgd.OutputFile(ab_output, extensions=['.yaml'])
            },
            mgd.TempSpace("combined")
        )
    )

    return workflow
