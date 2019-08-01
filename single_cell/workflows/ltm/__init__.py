'''
Created on July 31, 2018

@author: pwalters
'''

import pypeliner.managed as mgd

import pypeliner


def create_ltm_workflow(hmmcopy,
                        cn_matrix,
                        output_gml,
                        output_rooted_gml,
                        cnv_annots_csv,
                        cnv_tree_edges_csv,
                        cnv_data_csv,
                        output_rmd,
                        config,
                        root_id,
                        root_id_file,
                        number_jobs,
                        ploidy):
    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('timepoint'),
        value=list(hmmcopy.keys()),
    )

    workflow.transform(
        name='generate_cn_matrices',
        axes=('timepoint',),
        ctx={'mem': config['memory']['med'], 'pool_id': config['pools']['standard'], 'ncpus': 1},
        func='single_cell.workflows.ltm.tasks.generate_cn_matrices',
        args=(
            mgd.InputFile('hmmcopy.h5', 'timepoint', fnames=hmmcopy),
            mgd.TempOutputFile('cn_matrix.csv', 'timepoint'),
            str(ploidy),
        ),
    )

    # Generate copy number matrix
    workflow.transform(
        name='generate_cn_matrix',
        ctx={'mem': config['memory']['med'], 'pool_id': config['pools']['standard'], 'ncpus': 1},
        func='single_cell.workflows.ltm.tasks.merge_cn_matrices',
        args=(
            mgd.TempInputFile('cn_matrix.csv', 'timepoint'),
            mgd.OutputFile(cn_matrix),
        ),
    )

    node_pair_csvs = []
    for job in range(number_jobs):
        node_pair_csvs.append('list_{}.csv'.format(job))

    workflow.transform(
        name='generate_input_csvs',
        ctx={'mem': config['memory']['med'], 'pool_id': config['pools']['standard'], 'ncpus': 1},
        func='single_cell.workflows.ltm.tasks.generate_node_pair_csvs',
        args=(
            mgd.InputFile(cn_matrix),
            number_jobs,
            [mgd.TempOutputFile(csv) for csv in node_pair_csvs],
        ),
    )

    distance_csvs = []
    for job in range(number_jobs):
        distance_csvs.append('distance_list_{}.csv'.format(job))

    workflow.transform(
        name='calculate_distances',
        ctx={'mem': config['memory']['med'], 'pool_id': config['pools']['standard'], 'ncpus': 1},
        func='single_cell.workflows.ltm.tasks.calculate_distances',
        args=(
            [mgd.TempInputFile(csv) for csv in node_pair_csvs],
            mgd.InputFile(cn_matrix),
            [mgd.TempOutputFile(csv) for csv in distance_csvs],
            config,
        ),
    )

    # Generates a minimum spanning tree
    workflow.transform(
        name='generate_tree',
        ctx={'mem': config['memory']['med'], 'pool_id': config['pools']['standard'], 'ncpus': 1},
        func='single_cell.workflows.ltm.scripts.learn_CL_from_distance.learn_CL_from_distance',
        args=(
            [mgd.TempInputFile(csv) for csv in distance_csvs],
            mgd.OutputFile(output_gml),
        ),
    )

    workflow.transform(
        name='generate_cellscape_inputs',
        ctx={'mem': config['memory']['med'], 'pool_id': config['pools']['standard'], 'ncpus': 1},
        func='single_cell.workflows.ltm.tasks.generate_cellscape_inputs',
        args=(
            mgd.InputFile(cn_matrix),
            mgd.OutputFile(cnv_annots_csv),
            mgd.OutputFile(cnv_tree_edges_csv),
            mgd.OutputFile(cnv_data_csv),
            mgd.InputFile(output_gml),
            mgd.OutputFile(output_rooted_gml),
            root_id,
            mgd.OutputFile(root_id_file),
        ),
    )

    workflow.transform(
        name='create_cellscape_rmarkdown',
        ctx={'mem': config['memory']['med'], 'pool_id': config['pools']['standard'], 'ncpus': 1},
        func='single_cell.workflows.ltm.tasks.move_cellscape',
        args=(
            mgd.OutputFile(output_rmd),
        ),
    )

    return workflow
