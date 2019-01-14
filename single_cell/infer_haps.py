'''
Created on Aug 29, 2018

@author: dgrewal
'''

import os
import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers
from workflows import extract_seqdata
import single_cell

def infer_haps_workflow(workflow, args):

    config = helpers.load_config(args)
    config = config['infer_haps']

    baseimage = config['docker']['single_cell_pipeline']

    haps_dir = os.path.join(args["out_dir"], "infer_haps")

    haplotypes_filename = os.path.join(haps_dir, "results", "haplotypes.tsv")
    allele_counts_filename = os.path.join(haps_dir, "results", "allele_counts.tsv")

    workflow.setobj(
        obj=mgd.OutputChunks('chromosome'),
        value=config['chromosomes']
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
                    fnames=bam_files,
                    extensions=['.bai']
                ),
                mgd.TempOutputFile("tumour.h5", "cell_id"),
                config.get('extract_seqdata', {}),
                config['ref_data_dir'],
                config,
            )
        )

        workflow.transform(
            name='merge_all_seqdata',
            ctx={'mem_retry_increment': 2, 'ncpus': 1, 'docker_image': baseimage},
            func="single_cell.workflows.titan.tasks.merge_overlapping_seqdata",
            args=(
                mgd.TempOutputFile("seqdata_normal_all_cells_merged.h5"),
                mgd.TempInputFile("tumour.h5", "cell_id"),
                config["chromosomes"]
            ),
        )
    else:
        workflow.subworkflow(
            name="extract_seqdata",
            func=extract_seqdata.create_extract_seqdata_workflow,
            args=(
                mgd.InputFile(args['input_bam']),
                mgd.TempOutputFile("seqdata_normal_all_cells_merged.h5"),
                config.get('extract_seqdata', {}),
                config['ref_data_dir'],
                config,
            ),
            kwargs={'multiprocess': True}
        )

    if args["normal"]:
        workflow.transform(
            name='infer_snp_genotype',
            axes=('chromosome',),
            ctx={'mem': 16, 'mem_retry_increment': 2, 'ncpus': 1, 'docker_image': baseimage},
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
            ctx={'mem': 16, 'mem_retry_increment': 2, 'ncpus': 1, 'docker_image': baseimage},
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
        ctx={'mem': 16, 'mem_retry_increment': 2, 'ncpus': 1, 'docker_image': baseimage},
        func='remixt.analysis.haplotype.infer_haps',
        args=(
            mgd.TempOutputFile('haps.tsv', 'chromosome'),
            mgd.TempInputFile('snp_genotype.tsv', 'chromosome'),
            mgd.InputInstance('chromosome'),
            mgd.TempSpace('haplotyping', 'chromosome'),
            config,
            config['ref_data_dir'],
        )
    )

    workflow.transform(
        name='merge_haps',
        ctx={'mem': 16, 'mem_retry_increment': 2, 'ncpus': 1, 'docker_image': baseimage},
        func='remixt.utils.merge_tables',
        args=(
            mgd.OutputFile(haplotypes_filename),
            mgd.TempInputFile('haps.tsv', 'chromosome'),
        )
    )

    workflow.transform(
        name='create_segments',
        ctx={'mem': 20, 'mem_retry_increment': 2, 'ncpus': 1, 'docker_image': baseimage},
        func='remixt.analysis.segment.create_segments',
        args=(
            mgd.TempOutputFile('segments.tsv'),
            config,
            config['ref_data_dir'],
        ),
    )

    workflow.transform(
        name='haplotype_allele_readcount',
        ctx={'mem': 20, 'mem_retry_increment': 2, 'ncpus': 1, 'docker_image': baseimage},
        func='remixt.analysis.readcount.haplotype_allele_readcount',
        args=(
            mgd.OutputFile(allele_counts_filename),
            mgd.TempInputFile('segments.tsv'),
            mgd.TempInputFile('seqdata_normal_all_cells_merged.h5'),
            mgd.InputFile(haplotypes_filename),
            config,
        ),
    )

    info_file = os.path.join(args["out_dir"],'results','infer_haps', "info.yaml")

    results = {
        'infer_haps_allele_counts': helpers.format_file_yaml(allele_counts_filename),
        'infer_haps_data': helpers.format_file_yaml(haplotypes_filename),
    }

    if args['input_yaml']:
        input_datasets = {k: helpers.format_file_yaml(v) for k,v in bam_files.iteritems()}
    else:
        input_datasets = helpers.format_file_yaml(args['input_bam'])

    metadata = {
        'infer_haps': {
            'version': single_cell.__version__,
            'results': results,
            'containers': config['docker'],
            'input_datasets': input_datasets,
            'output_datasets': None
        }
    }

    workflow.transform(
        name='generate_meta_yaml',
        ctx={'mem': 4, 'mem_retry_increment': 2, 'ncpus': 1, 'docker_image': baseimage},
        func="single_cell.utils.helpers.write_to_yaml",
        args=(
            mgd.OutputFile(info_file),
            metadata
        )
    )

    return workflow


def infer_haps_from_bulk_normal(
    normal_bam_file,
    normal_seqdata_file,
    haplotypes_filename,
    config,
):
    remixt_config = config.get('extract_seqdata', {})
    remixt_ref_data_dir = config['ref_data_dir']

    chromosomes = config['chromosomes'] #remixt.config.get_chromosomes(config, remixt_ref_data_dir)

    workflow = pypeliner.workflow.Workflow()

    workflow.subworkflow(
        name='extract_seqdata',
        func='remixt.workflow.create_extract_seqdata_workflow',
        args=(
            mgd.InputFile(normal_bam_file, extensions=['.bai']),
            mgd.OutputFile(normal_seqdata_file),
            remixt_config,
            remixt_ref_data_dir,
        ),
    )

    workflow.setobj(
        obj=mgd.OutputChunks('chromosome'),
        value=chromosomes,
    )

    workflow.transform(
        name='infer_snp_genotype',
        axes=('chromosome',),
        ctx={'mem': 16},
        func='remixt.analysis.haplotype.infer_snp_genotype_from_normal',
        args=(
            mgd.TempOutputFile('snp_genotype.tsv', 'chromosome'),
            mgd.InputFile(normal_seqdata_file),
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
        ctx={'mem': 16},
        func='remixt.utils.merge_tables',
        args=(
            mgd.OutputFile(haplotypes_filename),
            mgd.TempInputFile('haplotypes.tsv', 'chromosome'),
        )
    )

    return workflow


def extract_allele_readcounts(
    haplotypes_filename,
    tumour_cell_bams,
    tumour_cell_seqdata,
    allele_counts_filename,
    config,
):
    tumour_cell_seqdata = {cell_id: tumour_cell_seqdata[cell_id] for cell_id in tumour_cell_bams}

    remixt_config = config.get('extract_seqdata', {})
    remixt_ref_data_dir = config['ref_data_dir']

    workflow = pypeliner.workflow.Workflow()

    workflow.set_filenames('cell.bam', 'cell_id', fnames=tumour_cell_bams)
    workflow.set_filenames('seqdata.h5', 'cell_id', fnames=tumour_cell_seqdata)

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=tumour_cell_bams.keys(),
    )

    workflow.subworkflow(
        name='create_chromosome_seqdata',
        axes=('cell_id',),
        func='remixt.workflow.create_extract_seqdata_workflow',
        args=(
            mgd.InputFile('cell.bam', 'cell_id'),
            mgd.OutputFile('seqdata.h5', 'cell_id', axes_origin=[]),
            remixt_config,
            remixt_ref_data_dir,
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
        ctx={'mem': 20},
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
        ctx={'mem': 20},
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

