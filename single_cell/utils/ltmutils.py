'''
Created on July 31, 2018

@author: pwalters
'''

import logging

import pandas as pd


def read_input_file(input_file):
    inputs = pd.read_csv(input_file, dtype=str)

    for column in ('timepoint', 'hmmcopy',):
        if column not in inputs.columns:
            raise Exception(
                'input_csv should contain {}'.format(column))

    timepoints = list(sorted(inputs['timepoint'].unique()))

    if inputs.duplicated(['timepoint']).any():
        raise Exception('duplicate timepoints in input_csv')

    hmmcopy = dict()
    for _, row in inputs.iterrows():
        hmmcopy[row['timepoint']] = row['hmmcopy'].strip()

    return hmmcopy, timepoints


def get_cn_matrix_from_hdf(hmmcopy_hdf_file, ploidy='0'):
    df = pd.read_hdf(hmmcopy_hdf_file, '/hmmcopy/reads/' + ploidy)

    df["bin"] = list(zip(df.chr, df.start, df.end))
    df = df.pivot(index='cell_id', columns='bin', values='state')
    chromosomes = map(str, range(1, 23)) + ['X', 'Y']
    bins = pd.DataFrame(df.columns.values.tolist(),
                        columns=['chr', 'start', 'end'])
    bins["chr"] = pd.Categorical(bins["chr"], chromosomes)
    bins = bins.sort_values(['start', ])
    bins = [tuple(v) for v in bins.values.tolist()]
    df = df.sort_values(bins, axis=0).T

    dropped_cells = df.columns[df.isna().all()].tolist()

    if len(dropped_cells) != 0:
        logging.getLogger("single_cell.helpers.ltmutils").warn(
            'Dropping {} cells: {}'.format(len(dropped_cells), dropped_cells)
        )

    df = df.loc[:, ~df.isna().all()].astype(int)
    df.columns = df.columns.astype(str)
    df = df.reset_index()

    chrom = []
    start = []
    end = []
    width = []
    for i, b in df['bin'].items():
        chrom.append(b[0])
        start.append(b[1])
        end.append(b[2])
        width.append(b[2] - b[1] + 1)
    df['chr'] = chrom
    df['start'] = start
    df['end'] = end
    df['width'] = width

    df = df.drop(columns='bin')

    return df, dropped_cells


def get_root(cells_list, root_id_file):
    for cell in cells_list:
        if 'SA928' in cell:
            with open(root_id_file, 'w') as outfile:
                outfile.write(cell + '\n')
            outfile.close()
            return cell

    raise Exception('No SA928 cells in the copy number matrix.')
