'''
Created on Apr 13, 2018

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
import biowrappers
import tasks

import remixt.workflow
from biowrappers.components.copy_number_calling import titan
import biowrappers.components.io.hdf5.tasks as hdf5_tasks

default_chromosomes = [str(a) for a in xrange(1, 23)] + ['X']



def create_extract_seqdata_workflow(
     bam_filename,
     seqdata_filename,
     config,
     ref_data_dir
):
    chromosomes = remixt.config.get_chromosomes(config, ref_data_dir)
    snp_positions_filename = remixt.config.get_filename(config, ref_data_dir, 'snp_positions')

    bam_max_fragment_length = remixt.config.get_param(config, 'bam_max_fragment_length')
    bam_max_soft_clipped = remixt.config.get_param(config, 'bam_max_soft_clipped')
    bam_check_proper_pair = remixt.config.get_param(config, 'bam_check_proper_pair')

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(obj=mgd.OutputChunks('chromosome'), value=chromosomes)

    workflow.transform(
        name='create_chromosome_seqdata',
        ctx={'mem': 16},
        func=tasks.create_chromosome_seqdata,
        args=(
            mgd.TempOutputFile('seqdata', 'chromosome', axes_origin=[]),
            mgd.InputFile(bam_filename),
            snp_positions_filename,
            chromosomes,
            bam_max_fragment_length,
            bam_max_soft_clipped,
            bam_check_proper_pair,
        ),
    )

    workflow.transform(
        name='merge_seqdata',
        ctx={'mem': 16},
        func=remixt.seqdataio.merge_seqdata,
        args=(
            mgd.OutputFile(seqdata_filename),
            mgd.TempInputFile('seqdata', 'chromosome'),
        ),
    )

    return workflow



def create_titan_workflow(
        normal_bam, normal_bai, tumour_bam, tumour_bai, ref_genome,
        out_file, raw_data_dir,
        config, args, tumour_cells, normal_cells):

    results_files = os.path.join(raw_data_dir, 'results', 'sample.h5')

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('tumour_cell_id'),
        value=tumour_cells,
    )

    workflow.setobj(
        obj=mgd.OutputChunks('normal_cell_id'),
        value=normal_cells,
    )


    workflow.subworkflow(
        name='extract_seqdata_workflow_normal',
        axes=('normal_cell_id',),
        func=create_extract_seqdata_workflow,
        args=(
            pypeliner.managed.InputFile('bam', 'normal_cell_id', fnames=normal_bam),
            pypeliner.managed.TempOutputFile('normal_sample.h5', 'normal_cell_id'),
            config['titan_params'].get('extract_seqdata', {}),
            config['titan_params']['ref_data_dir'],
        ),
    )

    workflow.subworkflow(
        name='extract_seqdata_workflow_tumour',
        axes=('tumour_cell_id',),
        func=create_extract_seqdata_workflow,
        args=(
            pypeliner.managed.InputFile('bam', 'tumour_cell_id', fnames=tumour_bam),
            pypeliner.managed.TempOutputFile('tumour_sample.h5', 'tumour_cell_id'),
            config['titan_params'].get('extract_seqdata', {}),
            config['titan_params']['ref_data_dir'],
        ),
    )

#     normal_seqdata = "/temp/titan_workflow/normal_cell_id/{normal_cell_id}/normal_sample.h5"
#     tumour_seqdata = "/temp/titan_workflow/tumour_cell_id/{tumour_cell_id}/tumour_sample.h5"

    workflow.transform(
        name='prepare_normal_data',
        ctx={'mem': 16, 'num_retry': 3, 'mem_retry_increment': 4},
        axes=('normal_cell_id',),
        func=titan.tasks.prepare_normal_data,
        args=(
#               pypeliner.managed.InputFile('normal_sample.h5', 'normal_cell_id', template=normal_seqdata),
            pypeliner.managed.TempInputFile('normal_sample.h5', 'normal_cell_id'),
            pypeliner.managed.TempOutputFile('normal.wig', 'normal_cell_id'),
            pypeliner.managed.TempOutputFile('het_positions.tsv', 'normal_cell_id'),
            config["titan_params"],
        ),
    )

    workflow.transform(
        name='merge_het_positions',
        ctx={'mem': 16, 'num_retry': 3, 'mem_retry_increment': 4},
        func=tasks.merge_het_positions,
        args=(
            pypeliner.managed.TempInputFile('het_positions.tsv', 'normal_cell_id'),
            pypeliner.managed.TempOutputFile('het_positions.tsv'),
        ),
    )


    workflow.transform(
        name='prepare_tumour_data',
        axes=('tumour_cell_id',),
        ctx={'mem': 20},
        func=titan.tasks.prepare_tumour_data,
        args=(
#               pypeliner.managed.InputFile('tumour_sample.h5', 'tumour_cell_id', template=tumour_seqdata),
            pypeliner.managed.TempInputFile('tumour_sample.h5', 'tumour_cell_id'),
            pypeliner.managed.TempInputFile('het_positions.tsv'),
            pypeliner.managed.TempOutputFile('tumour.wig', 'tumour_cell_id'),
            pypeliner.managed.TempOutputFile('tumour_alleles.tsv', 'tumour_cell_id'),
            config["titan_params"],
        ),
    )

    workflow.transform(
        name='merge_tumour_alleles',
        ctx={'mem': 16, 'num_retry': 3, 'mem_retry_increment': 4},
        func=tasks.merge_tumour_alleles,
        args=(
            pypeliner.managed.TempInputFile('tumour_alleles.tsv', 'tumour_cell_id'),
            pypeliner.managed.TempOutputFile('tumour_alleles.tsv'),
        ),
    )

    workflow.transform(
        name='merge_wigs_normal',
        ctx={'mem': 16, 'num_retry': 3, 'mem_retry_increment': 4},
        func=tasks.merge_wig_files,
        args=(
            pypeliner.managed.TempInputFile('normal.wig', 'normal_cell_id'),
            pypeliner.managed.TempOutputFile('normal.wig'),
        ),
    )

    workflow.transform(
        name='merge_wigs_tumour',
        ctx={'mem': 16, 'num_retry': 3, 'mem_retry_increment': 4},
        func=tasks.merge_wig_files,
        args=(
            pypeliner.managed.TempInputFile('tumour.wig', 'tumour_cell_id'),
            pypeliner.managed.TempOutputFile('tumour.wig'),
        ),
    )

    workflow.transform(
        name='create_intialization_parameters',
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
        func=titan.tasks.create_intialization_parameters,
        ret=pypeliner.managed.TempOutputObj('init_params', 'init_param_id'),
        args=(config["titan_params"],),
    )


    workflow.transform(
        name='run_titan',
        axes=('init_param_id',),
        ctx={'mem': 16, 'num_retry': 3, 'mem_retry_increment': 4},
        func=titan.tasks.run_titan,
        args=(
            pypeliner.managed.TempInputObj('init_params', 'init_param_id'),
            pypeliner.managed.TempInputFile('normal.wig'),
            pypeliner.managed.TempInputFile('tumour.wig'),
            pypeliner.managed.TempInputFile('tumour_alleles.tsv'),
            pypeliner.managed.TempOutputFile('cn.tsv', 'init_param_id'),
            pypeliner.managed.TempOutputFile('params.tsv', 'init_param_id'),
            config["titan_params"],
        ),
    )

    workflow.transform(
        name='select_solution',
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
        func=titan.tasks.select_solution,
        args=(
            pypeliner.managed.TempInputObj('init_params', 'init_param_id'),
            pypeliner.managed.TempInputFile('cn.tsv', 'init_param_id'),
            pypeliner.managed.TempInputFile('params.tsv', 'init_param_id'),
            pypeliner.managed.OutputFile('results', template=results_files),
            pypeliner.managed.OutputFile(os.path.join(raw_data_dir, 'output', 'cn_loci.tsv')),
            pypeliner.managed.OutputFile(
                os.path.join(raw_data_dir, 'output', 'cn_segments.tsv')),
            pypeliner.managed.OutputFile(os.path.join(raw_data_dir, 'output', 'cn_igv.tsv')),
            pypeliner.managed.OutputFile(os.path.join(raw_data_dir, 'output', 'params.tsv')),
            config,
            "MERGED"
        ),
        kwargs={
            'breakpoints_filename': None,
        },
    )

#     workflow.setobj(
#         obj=pypeliner.managed.OutputChunks('sample_id', 'chromosome'),
#         value=config.get('chromosomes', default_chromosomes),
#         axes=('sample_id',)
#     )
#
#     workflow.commandline(
#         name='plot_chromosome',
#         axes=('sample_id', 'chromosome'),
#         ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
#         args=(
#             'plot_titan_chromosome.R',
#             pypeliner.managed.Instance('chromosome'),
#             pypeliner.managed.InputFile(os.path.join(raw_data_dir, 'output', '{sample_id}_cn_loci.tsv'), 'sample_id'),
#             pypeliner.managed.InputFile(os.path.join(raw_data_dir, 'output', '{sample_id}_params.tsv'), 'sample_id'),
#             pypeliner.managed.OutputFile(
#                 os.path.join(raw_data_dir, 'output', '{sample_id}_chr_{chromosome}.png'), 'sample_id', 'chromosome'),
#         ),
#     )
#
#     workflow.transform(
#         name='merge_results',
#         ctx={'mem': 8, 'num_retry': 3, 'mem_retry_increment': 2},
#         func=hdf5_tasks.merge_hdf5,
#         args=(
#             pypeliner.managed.InputFile('results', 'sample_id', template=results_files),
#             pypeliner.managed.OutputFile(out_file),
#         ),
#         kwargs={
#             'table_names': '/sample_{}',
#         },
#     )



    return workflow
