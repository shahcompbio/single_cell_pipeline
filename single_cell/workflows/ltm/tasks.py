'''
Created on July 31, 2018

@author: pwalters
'''

import multiprocessing
import os

import pandas as pd
from scipy.special import comb
from single_cell.utils import helpers
from single_cell.utils import ltmutils

import pypeliner

scripts_directory = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)),
    'scripts')


def generate_cn_matrices(hmmcopy, cn_matrix, ploidy='0'):
    cn_df, dropped_cells = ltmutils.get_cn_matrix_from_hdf(hmmcopy, ploidy)
    metrics_summary_df = pd.read_hdf(hmmcopy, '/hmmcopy/metrics/' + ploidy)
    reads_df = pd.read_hdf(hmmcopy, '/hmmcopy/reads/' + ploidy)

    if dropped_cells:
        metrics_summary_df = metrics_summary_df.drop(
            index=metrics_summary_df.loc[metrics_summary_df['cell_id'].isin(dropped_cells)].index)

    # Filter out cells with quality less than 0.75
    cn_df = cn_df.drop(columns=metrics_summary_df[metrics_summary_df['quality'] < 0.75]['cell_id'])

    # Filter out bins with mappability < 0.99
    i = reads_df.loc[(reads_df['cell_id'] == reads_df.iloc[0]['cell_id']) & (reads_df['map'] >= 0.99)].index
    cn_df = cn_df.iloc[i]

    cn_df.to_csv(cn_matrix, index=False)


def merge_cn_matrices(infiles, outfile):
    cn_list = []

    for cn_matrix in infiles.values():
        cn_list.append(pd.read_csv(cn_matrix))

    cn_matrix = pd.concat(cn_list, axis=1, join='outer')
    cn_matrix = cn_matrix.T.drop_duplicates().T  # Get rid of any duplicate columns
    cn_matrix = cn_matrix.set_index(['chr', 'start', 'end', 'width'])

    cn_matrix.to_csv(outfile)


def generate_node_pair_csvs(cn_matrix, desired_job_no, node_pair_csvs):
    cells_list = pd.read_csv(cn_matrix).columns.tolist()[4:]

    total_edges = comb(len(cells_list), 2)
    edg_per_node = int(total_edges / float(desired_job_no))

    tmp_counter = 0
    global_counter = 0
    file_counter = 1

    outfile = open(node_pair_csvs[0], 'w')

    for i in range(len(cells_list)):
        for j in range(i + 1, len(cells_list)):
            outfile.write(cells_list[i] + ',' + cells_list[j] + '\n')

            tmp_counter += 1
            global_counter += 1

            if tmp_counter == edg_per_node and global_counter < total_edges and file_counter < desired_job_no:
                tmp_counter = 0
                outfile.close()

                outfile = open(node_pair_csvs[file_counter], 'w')
                file_counter += 1

            if global_counter == total_edges:
                outfile.close()


def _calculate_distances_worker(node_pair_csv, outfile, cn_matrix):
    script = os.path.join(scripts_directory, 'calculate_distance.py')

    cmd = ['python', script, '-input_file', node_pair_csv,
           '-output_path', outfile, '-data_path', cn_matrix]

    helpers.run_cmd(cmd)


def calculate_distances(node_pair_csvs, cn_matrix, outfiles, config):
    count = config.get('threads', multiprocessing.cpu_count())
    pool = multiprocessing.Pool(processes=count)

    tasks = []

    for job in zip(node_pair_csvs, outfiles):
        node_pair_csv = job[0]
        outfile = job[1]

        task = pool.apply_async(_calculate_distances_worker,
                                args=(node_pair_csv, outfile, cn_matrix))

        tasks.append(task)

    pool.close()
    pool.join()

    [task.get() for task in tasks]


## VISUALIZATION ##

def generate_cellscape_inputs(cn_matrix, annotations, edges_list, cn_data, tree_gml, rooted_tree_gml, root_id,
                              root_id_file):
    cells_list = pd.read_csv(cn_matrix).columns.tolist()[4:]

    if root_id:
        if root_id not in cells_list:
            raise Exception('Root ID {root_id} is not a cell in the copy number matrix.'.format(root_id=root_id))
    else:
        root_id = ltmutils.get_root(cells_list, root_id_file)

    script = os.path.join(scripts_directory, 'generate_cellscape_inputs.py')

    cmd = ['python', script,
           '--path_to_data', cn_matrix,
           '--path_to_annotations', annotations,
           '--path_to_edges_list', edges_list,
           '--path_to_cn_data', cn_data,
           '--path_to_tree', tree_gml,
           '--path_to_rooted_tree', rooted_tree_gml,
           '--root_id', root_id]

    pypeliner.commandline.execute(*cmd)


def move_cellscape(outfile):
    cellscape_rmarkdown = os.path.join(scripts_directory, 'cellscape.Rmd')
    helpers.copy_file(cellscape_rmarkdown, outfile)
