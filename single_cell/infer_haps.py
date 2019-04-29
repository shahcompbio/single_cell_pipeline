'''
Created on Aug 29, 2018

@author: dgrewal
'''

import os
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers


def infer_haps(
    bam_file,
    seqdata_file,
    haplotypes_filename,
    allele_counts_filename,
    config,
    normal=False,
):

    baseimage = {'docker_image': config['docker']['single_cell_pipeline']}

    remixt_config = config.get('extract_seqdata', {})
    remixt_ref_data_dir = config['ref_data_dir']

    chromosomes = config['chromosomes']

    ctx = dict(mem_retry_increment=2, disk_retry_increment=50, ncpus=1, **baseimage)
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    if isinstance(bam_file, dict):
        workflow.setobj(
            obj=mgd.OutputChunks('cell_id'),
            value=list(bam_file.keys()),
        )

        # dont parallelize over chromosomes for per cell bams
        workflow.subworkflow(
            name="extract_seqdata",
            axes=('cell_id',),
            func='single_cell.workflows.extract_seqdata.create_extract_seqdata_workflow',
            args=(
                mgd.InputFile(
                    'bam_markdups',
                    'cell_id',
                    fnames=bam_file,
                    extensions=['.bai']
                ),
                mgd.TempOutputFile('seqdata_cell.h5','cell_id'),
                config.get('extract_seqdata', {}),
                config['ref_data_dir'],
                config,
            )
        )
        workflow.transform(
            name='merge_all_seqdata',
            func="single_cell.workflows.titan.tasks.merge_overlapping_seqdata",
            args=(
                mgd.OutputFile(seqdata_file),
                mgd.TempInputFile("seqdata_cell.h5", "cell_id"),
                config["chromosomes"]
            ),
        )

    else:
        # if its a single bam, then its probably whole genome
        # so parallelize over chromosomes
        workflow.subworkflow(
            name='extract_seqdata',
            func='remixt.workflow.create_extract_seqdata_workflow',
            ctx={'disk': 150},
            args=(
                mgd.InputFile(bam_file, extensions=['.bai']),
                mgd.OutputFile(seqdata_file),
                remixt_config,
                remixt_ref_data_dir,
            ),
        )

    workflow.setobj(
        obj=mgd.OutputChunks('chromosome'),
        value=chromosomes,
    )

    if normal:
        func = 'remixt.analysis.haplotype.infer_snp_genotype_from_normal'
    else:
        func = 'remixt.analysis.haplotype.infer_snp_genotype_from_tumour'

    workflow.transform(
        name='infer_snp_genotype',
        axes=('chromosome',),
        ctx=dict(mem=16, **ctx),
        func=func,
        args=(
            mgd.TempOutputFile('snp_genotype.tsv', 'chromosome'),
            mgd.InputFile(seqdata_file),
            mgd.InputInstance('chromosome'),
            config,
        ),
    )

    workflow.transform(
        name='infer_haps',
        axes=('chromosome',),
        ctx=dict(mem=16, **ctx),
        func='remixt.analysis.haplotype.infer_haps',
        args=(
            mgd.TempOutputFile('haplotypes.tsv', 'chromosome'),
            mgd.TempInputFile('snp_genotype.tsv', 'chromosome'),
            mgd.InputInstance('chromosome'),
            mgd.TempSpace('haplotyping', 'chromosome'),
            remixt_config,
            remixt_ref_data_dir,
        ),
    )

    workflow.transform(
        name='merge_haps',
        ctx=dict(mem=16, **ctx),
        func='remixt.utils.merge_tables',
        args=(
            mgd.OutputFile(haplotypes_filename),
            mgd.TempInputFile('haplotypes.tsv', 'chromosome'),
        )
    )

    workflow.transform(
        name='create_segments',
        ctx=dict(mem=16, **ctx),
        func='remixt.analysis.segment.create_segments',
        args=(
            mgd.TempOutputFile('segments.tsv'),
            config,
            config['ref_data_dir'],
        ),
    )

    workflow.transform(
        name='haplotype_allele_readcount',
        ctx=dict(mem=16, **ctx),
        func='remixt.analysis.readcount.haplotype_allele_readcount',
        args=(
            mgd.OutputFile(allele_counts_filename),
            mgd.TempInputFile('segments.tsv'),
            mgd.InputFile(seqdata_file),
            mgd.InputFile(haplotypes_filename),
            config,
        ),
    )

    return workflow


def extract_allele_readcounts(
    haplotypes_filename,
    cell_bams,
    cell_seqdata,
    allele_counts_filename,
    config,
):
    baseimage = {'docker_image': config['docker']['single_cell_pipeline']}

    cell_seqdata = {cell_id: cell_seqdata[cell_id] for cell_id in cell_bams}

    remixt_config = config.get('extract_seqdata', {})
    remixt_ref_data_dir = config['ref_data_dir']

    workflow = pypeliner.workflow.Workflow(ctx=baseimage)

    workflow.set_filenames('cell.bam', 'cell_id', fnames=cell_bams)
    workflow.set_filenames('seqdata.h5', 'cell_id', fnames=cell_seqdata)

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=list(cell_bams.keys()),
    )

    workflow.subworkflow(
        name='create_chromosome_seqdata',
        axes=('cell_id',),
        func='single_cell.workflows.extract_seqdata.create_extract_seqdata_workflow',
        args=(
            mgd.InputFile('cell.bam', 'cell_id', extensions=['.bai']),
            mgd.OutputFile('seqdata.h5', 'cell_id', axes_origin=[]),
            config.get('extract_seqdata', {}),
            config['ref_data_dir'],
            config,
        ),
    )

    #TODO Segments with bin width from single cell
    workflow.transform(
        name='create_segments',
        func='remixt.analysis.segment.create_segments',
        args=(
            mgd.TempOutputFile('segments.tsv'),
            remixt_config,
            remixt_ref_data_dir,
        ),
    )

    workflow.transform(
        name='haplotype_allele_readcount',
        axes=('cell_id',),
        ctx={'mem': 16},
        func='remixt.analysis.readcount.haplotype_allele_readcount',
        args=(
            mgd.TempOutputFile('allele_counts.tsv', 'cell_id', axes_origin=[]),
            mgd.TempInputFile('segments.tsv'),
            mgd.InputFile('seqdata.h5', 'cell_id'),
            mgd.InputFile(haplotypes_filename),
            remixt_config,
        ),
    )

    workflow.transform(
        name='merge_allele_readcount',
        ctx={'mem': 16},
        func='single_cell.utils.csvutils.concatenate_csv',
        args=(
            mgd.TempInputFile('allele_counts.tsv', 'cell_id'),
            mgd.OutputFile(allele_counts_filename),
        ),
        kwargs={
            'key_column': 'cell_id',
            'sep': '\t',
        },
    )

    return workflow



def infer_haps_workflow(args):
    config = helpers.load_config(args)
    config = config['infer_haps']
    baseimage = config['docker']['single_cell_pipeline']

    ctx = dict(mem_retry_increment=2, disk_retry_increment=50, ncpus=1, baseimage=baseimage)
    workflow = pypeliner.workflow.Workflow(ctx=ctx)


    haps_dir = os.path.join(args["out_dir"], "infer_haps")
    seqdata_merged = os.path.join(haps_dir, "results", "seqdata.h5")
    haplotypes_filename = os.path.join(haps_dir, "results", "haplotypes.tsv")
    allele_counts_filename = os.path.join(haps_dir, "results", "allele_counts.tsv")

    data = helpers.load_pseudowgs_input(args['input_yaml'])
    tumour_wgs = data['tumour_wgs']
    normal_wgs = data['normal_wgs']
    tumour_cells = data['tumour_cells']
    normal_cells = data['normal_cells']

    if args['normal']:
        bam_file = normal_cells if normal_cells else normal_wgs
    else:
        bam_file = tumour_cells if tumour_cells else tumour_wgs

    if isinstance(bam_file, dict):
        workflow.setobj(
            obj=mgd.OutputChunks('cell_id'),
            value=list(bam_file.keys()),
        )
        bam_file = mgd.InputFile('tumour.bam', 'cell_id',fnames=bam_file, extensions=['.bai'])
    else:
        bam_file = mgd.InputFile(bam_file, extensions=['.bai'])


    workflow.subworkflow(
        name='infer_haps',
        func=infer_haps,
        args=(
            bam_file,
            mgd.OutputFile(seqdata_merged),
            mgd.OutputFile(haplotypes_filename),
            mgd.OutputFile(allele_counts_filename),
            config,
        ),
        kwargs={'normal': args['normal']},
    )

    return workflow
