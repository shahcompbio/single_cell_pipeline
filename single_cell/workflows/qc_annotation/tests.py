import os

import pandas as pd
from single_cell.utils import csvutils
from single_cell.workflows.qc_annotation import tasks


def test_contamination(tmpdir):
    data = {}

    cols = [
        'fastqscreen_nohit',
        'fastqscreen_grch37',
        'fastqscreen_grch37_multihit',
        'fastqscreen_mm10',
        'fastqscreen_mm10_multihit',
        'fastqscreen_salmon',
        'fastqscreen_salmon_multihit'
    ]

    for i in range(5):
        data[i] = {'cell_id': 'SA123_A123_R{0}_C{0}'.format(i)}
        for col in cols:
            data[i][col] = i * 10
        data[i]['fastqscreen_grch37'] = i * 1000
        data[i]['fastqscreen_mm10'] = i * 100

    for i in range(5, 10):
        data[i] = {'cell_id': 'SA123_A123_R{0}_C{0}'.format(i)}
        for col in cols:
            data[i][col] = (i * 10)
        data[i]['fastqscreen_grch37'] = i * 1000

    data = pd.DataFrame.from_dict(data, orient='index')
    data['total_reads'] = data[cols].sum(axis=1)

    dtypes = {col: 'int' for col in cols}
    dtypes['cell_id'] = 'str'
    dtypes['total_reads'] = 'int'

    infile = os.path.join(tmpdir, 'input.csv.gz')
    outfile = os.path.join(tmpdir, 'output.csv.gz')

    csvutils.write_dataframe_to_csv_and_yaml(data, infile, dtypes)

    config = {'genomes': [{'name': 'grch37'}, {'name': 'mm10'}, {'name': 'salmon'}]}

    tasks.add_contamination_status(infile, outfile, config)

    output = csvutils.read_csv_and_yaml(outfile)

    assert output['is_contaminated'].tolist() == [False] + [True] * 4 + [False] * 5
