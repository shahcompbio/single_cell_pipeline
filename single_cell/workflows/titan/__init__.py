'''
Created on Apr 13, 2018

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd


default_chromosomes = [str(a) for a in xrange(1, 23)] + ['X']


def create_titan_workflow(normal_seqdata, tumour_seqdata, ref_genome,
                          raw_data_dir, out_file, config, args,
                          tumour_cells, normal_cells, cloneid):

    results_files = os.path.join(raw_data_dir, 'results', 'sample.h5')
    tumour_alleles_file = os.path.join(raw_data_dir, 'results', 'het_counts.h5')

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('tumour_cell_id'),
        value=tumour_cells,
    )

    workflow.setobj(
        obj=mgd.OutputChunks('normal_cell_id'),
        value=normal_cells,
    )

    workflow.transform(
        name='merge_all_normal_seqdata',
        ctx={'mem': config["memory"]['high'], 'pool_id': config['pools']['highmem'], 'ncpus':1},
        func="single_cell.workflows.titan.tasks.merge_overlapping_seqdata",
        args=(
            mgd.TempOutputFile("seqdata_normal_all_cells_merged.h5"),
            pypeliner.managed.InputFile(
                'normal_sample.h5',
                'normal_cell_id',
                fnames=normal_seqdata),
        ),
    )


    workflow.transform(
        name='prepare_normal_data',
        ctx={'mem': config["memory"]['high'],
             'pool_id': config['pools']['highmem'],
             'ncpus':1, 'num_retry': 3,
             'mem_retry_increment': 4},
        func="biowrappers.components.copy_number_calling.titan.tasks.prepare_normal_data",
        args=(
            mgd.TempInputFile("seqdata_normal_all_cells_merged.h5"),
            pypeliner.managed.TempOutputFile('normal.wig'),
            pypeliner.managed.TempOutputFile('het_positions.tsv'),
            config["titan_params"],
        ),
    )

    workflow.transform(
        name='prepare_tumour_data',
        axes=('tumour_cell_id',),
        ctx={'mem': config["memory"]['high'],
             'pool_id': config['pools']['highmem'],
             'ncpus':1},
        func="biowrappers.components.copy_number_calling.titan.tasks.prepare_tumour_data",
        args=(
            pypeliner.managed.InputFile(
                'tumour_sample.h5',
                'tumour_cell_id',
                fnames=tumour_seqdata),
            pypeliner.managed.TempInputFile('het_positions.tsv'),
            pypeliner.managed.TempOutputFile('tumour.wig', 'tumour_cell_id'),
            pypeliner.managed.TempOutputFile(
                'tumour_alleles.tsv',
                'tumour_cell_id'),
            config["titan_params"],
        ),
    )

    workflow.transform(
        name='merge_tumour_alleles',
        ctx={'mem': config["memory"]['high'],
             'pool_id': config['pools']['highmem'],
             'ncpus':1, 'num_retry': 3,
             'mem_retry_increment': 4},
        func="single_cell.workflows.titan.tasks.merge_tumour_alleles",
        args=(
            pypeliner.managed.TempInputFile(
                'tumour_alleles.tsv',
                'tumour_cell_id'),
            pypeliner.managed.TempOutputFile('tumour_alleles.tsv'),
        ),
    )

    workflow.transform(
        name='concat_tumour_alleles',
        ctx={'mem': config["memory"]['high'],
             'pool_id': config['pools']['highmem'],
             'ncpus':1, 'num_retry': 3,
             'mem_retry_increment': 4},
        func="single_cell.workflows.titan.tasks.concat_tumour_alleles",
        args=(
            pypeliner.managed.TempInputFile(
                'tumour_alleles.tsv',
                'tumour_cell_id'),
            pypeliner.managed.OutputFile(tumour_alleles_file),
        ),
    )

    workflow.transform(
        name='merge_wigs_tumour',
        ctx={'mem': config["memory"]['high'],
             'pool_id': config['pools']['highmem'],
             'ncpus':1, 'num_retry': 3,
             'mem_retry_increment': 4},
        func="single_cell.workflows.titan.tasks.merge_wig_files",
        args=(
            pypeliner.managed.TempInputFile('tumour.wig', 'tumour_cell_id'),
            pypeliner.managed.TempOutputFile('tumour.wig'),
        ),
    )

    workflow.transform(
        name='create_intialization_parameters',
        ctx={'mem': config["memory"]['low'],
             'pool_id': config['pools']['standard'],
             'ncpus':1, 'num_retry': 3,
             'mem_retry_increment': 2},
        func="biowrappers.components.copy_number_calling.titan.tasks.create_intialization_parameters",
        ret=pypeliner.managed.TempOutputObj('init_params', 'init_param_id'),
        args=(config["titan_params"],),
    )

    workflow.transform(
        name='run_titan',
        axes=('init_param_id',),
        ctx={'mem': config["memory"]['high'],
             'pool_id': config['pools']['highmem'],
             'ncpus':1, 'num_retry': 3,
             'mem_retry_increment': 4},
        func="biowrappers.components.copy_number_calling.titan.tasks.run_titan",
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
        ctx={'mem': config["memory"]['low'],
             'pool_id': config['pools']['standard'],
             'ncpus':1, 'num_retry': 3,
             'mem_retry_increment': 2},
        func="biowrappers.components.copy_number_calling.titan.tasks.select_solution",
        args=(
            pypeliner.managed.TempInputObj('init_params', 'init_param_id'),
            pypeliner.managed.TempInputFile('cn.tsv', 'init_param_id'),
            pypeliner.managed.TempInputFile('params.tsv', 'init_param_id'),
            pypeliner.managed.OutputFile('results', template=results_files),
            pypeliner.managed.OutputFile(
                os.path.join(
                    raw_data_dir,
                    'output',
                    'cn_loci.tsv')),
            pypeliner.managed.OutputFile(
                os.path.join(raw_data_dir, 'output', 'cn_segments.tsv')),
            pypeliner.managed.OutputFile(
                os.path.join(
                    raw_data_dir,
                    'output',
                    'cn_igv.tsv')),
            pypeliner.managed.OutputFile(
                os.path.join(
                    raw_data_dir,
                    'output',
                    'params.tsv')),
            config,
            cloneid
        ),
        kwargs={
            'breakpoints_filename': None,
        },
    )

    workflow.setobj(
        obj=mgd.OutputChunks('chromosome'),
        value=default_chromosomes,
    )

    workflow.commandline(
        name='plot_chromosome',
        axes=('chromosome',),
        ctx={'mem': config["memory"]['low'],
             'pool_id': config['pools']['standard'],
             'ncpus':1, 'num_retry': 3,
             'mem_retry_increment': 2},
        args=(
            'plot_titan_chromosome.R',
            pypeliner.managed.Instance('chromosome'),
            pypeliner.managed.InputFile(
                os.path.join(
                    raw_data_dir,
                    'output',
                    'cn_loci.tsv')),
            pypeliner.managed.InputFile(
                os.path.join(
                    raw_data_dir,
                    'output',
                    'params.tsv')),
            pypeliner.managed.OutputFile(
                os.path.join(raw_data_dir, 'output', 'chr_{chromosome}.png'), 'chromosome'),
        ),
    )

    #just leaving it here in case we parallelize by samples later.
    workflow.transform(
        name='merge_results',
        ctx={'mem': config["memory"]['low'],
             'pool_id': config['pools']['standard'],
             'ncpus':1, 'num_retry': 3,
             'mem_retry_increment': 2},
        func="biowrappers.components.io.hdf5.tasks.merge_hdf5",
        args=(
            {cloneid: pypeliner.managed.InputFile(
                'results', template=results_files)},
            pypeliner.managed.OutputFile(out_file),
        ),
        kwargs={
            'table_names': '/sample_{}'.format(cloneid),
        },
    )

    return workflow
