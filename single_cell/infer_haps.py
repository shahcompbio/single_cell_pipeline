'''
Created on Aug 29, 2018

@author: dgrewal
'''

import os
import pypeliner.managed as mgd
from single_cell.utils import helpers
from workflows import extract_seqdata


def infer_haps_workflow(workflow, args):

    config = helpers.load_config(args)
    remixt_config = config['titan_params'].get('extract_seqdata', {})

    singlecellimage = config['docker']['images']['single_cell_pipeline']
    ctx = {
              'mem_retry_increment': 2,
              'ncpus': 1,
              'image': singlecellimage['image'],
              'dockerize': config['docker']['dockerize'],
              'mounts': config['docker']['mounts'],
              'username': singlecellimage['username'],
              'password': singlecellimage['password'],
              'server': singlecellimage['server'],
          }

    haps_dir = os.path.join(args["out_dir"], "infer_haps")

    haplotypes_filename = os.path.join(haps_dir, "results", "haplotypes.tsv")
    allele_counts_filename = os.path.join(haps_dir, "results", "allele_counts.tsv")

    snp_positions_filename = remixt.config.get_filename(config, ref_data_dir, 'snp_positions')
    bam_max_fragment_length = remixt.config.get_param(config, 'bam_max_fragment_length')
    bam_max_soft_clipped = remixt.config.get_param(config, 'bam_max_soft_clipped')
    bam_check_proper_pair = remixt.config.get_param(config, 'bam_check_proper_pair')

    workflow.setobj(
        obj=mgd.OutputChunks('chromosome'),
        value=config['titan_params']['chromosomes']
    )

    if args['input_yaml']:
        bam_files, bai_files = helpers.get_bams(args['input_yaml'])
        cellids = helpers.get_samples(args['input_yaml'])

        workflow.setobj(
            obj=mgd.OutputChunks('cell_id'),
            value=cellids,
        )

        workflow.subworkflow(
            name="extract_seqdata",
            axes=('cell_id',),
            func=extract_seqdata.create_extract_seqdata_workflow,
            args=(
                mgd.InputFile(
                    'bam_markdups',
                    'cell_id',
                    fnames=bam_files),
                mgd.InputFile(
                    'bam_markdups_index',
                    'cell_id',
                    fnames=bai_files),
                mgd.TempOutputFile("tumour.h5", "cell_id"),
                config,
                config['titan_params'].get('extract_seqdata', {}),
                config['titan_params']['ref_data_dir'],
                snp_positions_filename,
                bam_max_fragment_length,
                bam_max_soft_clipped,
                bam_check_proper_pair,
            )
        )

        workflow.transform(
            name='merge_all_seqdata',
            ctx=dict(mem=config["memory"]['high'], pool_id=config['pools']['highmem'], **ctx),
            func="single_cell.workflows.titan.tasks.merge_overlapping_seqdata",
            args=(
                mgd.TempOutputFile("seqdata_normal_all_cells_merged.h5"),
                mgd.TempInputFile("tumour.h5", "cell_id"),
                config["titan_params"]["chromosomes"]
            ),
        )
    else:
        workflow.subworkflow(
            name="extract_seqdata",
            func=extract_seqdata.create_extract_seqdata_workflow,
            args=(
                mgd.InputFile(args['input_bam']),
                mgd.InputFile(args['input_bam'] + '.bai'),
                mgd.TempOutputFile("seqdata_normal_all_cells_merged.h5"),
                config,
                config['titan_params'].get('extract_seqdata', {}),
                config['titan_params']['ref_data_dir'],
                snp_positions_filename,
                bam_max_fragment_length,
                bam_max_soft_clipped,
                bam_check_proper_pair,
            ),
            kwargs={'multiprocess': True}
        )

    if args["normal"]:
        workflow.transform(
            name='infer_snp_genotype',
            axes=('chromosome',),
            ctx={'mem': 16},
            func='remixt.analysis.haplotype.infer_snp_genotype_from_normal',
            args=(
                mgd.TempOutputFile('snp_genotype.tsv', 'chromosome'),
                mgd.TempInputFile("seqdata_normal_all_cells_merged.h5"),
                mgd.InputInstance('chromosome'),
                config,
            ),
        )
    else:
        workflow.transform(
            name='infer_snp_genotype',
            axes=('chromosome',),
            ctx={'mem': 16},
            func='remixt.analysis.haplotype.infer_snp_genotype_from_tumour',
            args=(
                mgd.TempOutputFile('snp_genotype.tsv', 'chromosome'),
                {'sample': mgd.TempInputFile("seqdata_normal_all_cells_merged.h5")},
                mgd.InputInstance('chromosome'),
                config,
            ),
        )

    workflow.transform(
        name='infer_haps',
        axes=('chromosome',),
        ctx={'mem': 16},
        func='remixt.analysis.haplotype.infer_haps',
        args=(
            mgd.TempOutputFile('haps.tsv', 'chromosome'),
            mgd.TempInputFile('snp_genotype.tsv', 'chromosome'),
            mgd.InputInstance('chromosome'),
            mgd.TempSpace('haplotyping', 'chromosome'),
            config,
            config['titan_params']['ref_data_dir'],
        )
    )

    workflow.transform(
        name='merge_haps',
        ctx={'mem': 16},
        func='remixt.utils.merge_tables',
        args=(
            mgd.OutputFile(haplotypes_filename),
            mgd.TempInputFile('haps.tsv', 'chromosome'),
        )
    )

    workflow.transform(
        name='create_segments',
        func='remixt.analysis.segment.create_segments',
        args=(
            mgd.TempOutputFile('segments.tsv'),
            config,
            config['titan_params']['ref_data_dir'],
        ),
    )

    workflow.transform(
        name='haplotype_allele_readcount',
        ctx={'mem': 20},
        func='remixt.analysis.readcount.haplotype_allele_readcount',
        args=(
            mgd.OutputFile(allele_counts_filename),
            mgd.TempInputFile('segments.tsv'),
            mgd.TempInputFile('tumour.h5', 'cell_id'),
            mgd.InputFile(haplotypes_filename),
            config,
        ),
    )

    return workflow
