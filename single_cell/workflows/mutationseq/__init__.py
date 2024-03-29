'''
Created on Jul 24, 2017

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd
from single_cell.workflows.mutationseq.dtypes import dtypes


def museq_callback(record):
    return record.INFO['PR']


def create_museq_workflow(
        normal_bam, tumour_bam, snv_vcf, museq_csv,
        config):
    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1, 'num_retry': 3}

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('region'),
        value=list(normal_bam.keys()),
    )

    workflow.transform(
        name='subsample_tumour',
        ctx=dict(mem=config["memory"]['med']),
        axes=('region',),
        func='single_cell.workflows.mutationseq.tasks.subsample',
        args=(
            mgd.InputFile('merged_bam', 'region', fnames=tumour_bam, extensions=['.bai']),
            mgd.TempOutputFile("tumour_subsample.bam", "region", extensions=['.bai'])
        ),
    )

    workflow.transform(
        name='subsample_normal',
        ctx=dict(mem=config["memory"]['med']),
        axes=('region',),
        func='single_cell.workflows.mutationseq.tasks.subsample',
        args=(
            mgd.InputFile('normal.split.bam', 'region', fnames=normal_bam, extensions=['.bai']),
            mgd.TempOutputFile("normal_subsample.bam", "region", extensions=['.bai'])
        ),
    )

    workflow.transform(
        name='run_museq',
        ctx=dict(mem=config["memory"]['med']),
        axes=('region',),
        func='single_cell.workflows.mutationseq.tasks.run_museq',
        args=(
            mgd.TempInputFile("tumour_subsample.bam", "region", extensions=['.bai']),
            mgd.TempInputFile("normal_subsample.bam", "region", extensions=['.bai']),
            mgd.TempOutputFile('museq.vcf', 'region'),
            mgd.TempOutputFile('museq.log', 'region'),
            mgd.InputInstance('region'),
            config,
        ),
    )

    workflow.transform(
        name='finalise_region_vcfs',
        axes=('region',),
        ctx=dict(mem=config["memory"]['med']),
        func='biowrappers.components.io.vcf.tasks.finalise_vcf',
        args=(
            mgd.TempInputFile('museq.vcf', 'region'),
            mgd.TempOutputFile('museq.vcf.gz', 'region', extensions=['.tbi', '.csi']),
        ),
    )

    workflow.transform(
        name='merge_snvs',
        ctx=dict(mem=config["memory"]['med']),
        func='biowrappers.components.io.vcf.tasks.concatenate_vcf',
        args=(
            mgd.TempInputFile('museq.vcf.gz', 'region', extensions=['.tbi', '.csi']),
            mgd.TempOutputFile('museq.vcf.gz', extensions=['.tbi', '.csi']),
        ),
        kwargs={
            'allow_overlap': True,
        },
    )

    workflow.transform(
        name='finalise_vcf',
        func='biowrappers.components.io.vcf.tasks.finalise_vcf',
        ctx=dict(mem=config["memory"]['med']),
        args=(
            mgd.TempInputFile('museq.vcf.gz', extensions=['.tbi', '.csi']),
            mgd.OutputFile(snv_vcf, extensions=['.tbi', '.csi']),
        ),
    )

    workflow.transform(
        name='convert_museq_to_csv',
        func="biowrappers.components.io.vcf.tasks.convert_vcf_to_csv",
        ctx=ctx,
        args=(
            mgd.InputFile(snv_vcf),
            mgd.TempOutputFile('museq.csv'),
        ),
        kwargs={
            'score_callback': museq_callback,
        }
    )

    workflow.transform(
        name='prep_museq_csv',
        func='single_cell.utils.csvutils.rewrite_csv_file',
        args=(
            mgd.TempInputFile('museq.csv'),
            mgd.OutputFile(museq_csv, extensions=['.yaml'])
        ),
        kwargs={'dtypes': dtypes()['snv_museq']}
    )

    return workflow
