import re

import pandas as pd
from single_cell.utils import helpers


def group_lumpy_data(infile):
    key = None
    data = []

    with open(infile) as f:
        for line in f:
            if not line.startswith("\t"):
                if key is not None:
                    yield key, data
                key = tuple(line.rstrip().split())
                data = []
            else:
                data.append(line)
        if key is not None:
            yield key, data


def generate_primary_table(parsed_data):
    data = []

    for brk_call, brk_data in parsed_data:
        row_data = dict()

        row_data['chrom1'] = brk_call[0]
        row_data['start1'] = brk_call[1]
        row_data['end1'] = brk_call[2]

        row_data['chrom2'] = brk_call[3]
        row_data['start2'] = brk_call[4]
        row_data['end2'] = brk_call[5]

        row_data['breakpoint_id'] = brk_call[6]
        row_data['score'] = float(brk_call[7])
        row_data['strand1'] = brk_call[8]
        row_data['strand2'] = brk_call[9]
        row_data['type'] = brk_call[10].replace('TYPE:', '')

        brk_ids = re.split('[:;]', brk_call[11])
        for brk_id in brk_ids[1:]:
            sample, count = brk_id.split(',')
            row_data[sample] = int(count)

        row_data['strands'] = re.split('[:,]', brk_call[12])[1]

        max_brk = re.split('[:;]', brk_call[13])
        assert len(max_brk) == 5
        row_data['max_chr1'] = max_brk[1]
        row_data['max_pos1'] = max_brk[2]
        row_data['max_chr2'] = max_brk[3]
        row_data['max_pos2'] = max_brk[4]

        conf_interval = re.split('[:\-;]', brk_call[14])
        assert conf_interval[0] == '95'
        assert len(conf_interval) == 7
        row_data['confidence_interval_chr1'] = conf_interval[1]
        row_data['confidence_interval_start1'] = conf_interval[2]
        row_data['confidence_interval_end1'] = conf_interval[3]

        row_data['confidence_interval_chr2'] = conf_interval[4]
        row_data['confidence_interval_start2'] = conf_interval[5]
        row_data['confidence_interval_end2'] = conf_interval[6]

        data.append(row_data)

    df = pd.DataFrame(data)

    if df.empty:
        return df

    columns = df.columns.values

    order = [
        'breakpoint_id', 'chrom1', 'start1', 'end1', 'strand1',
        'max_chr1', 'max_pos1', 'confidence_interval_chr1',
        'confidence_interval_start1', 'confidence_interval_end1',
        'chrom2', 'start2', 'end2', 'strand2', 'max_chr2', 'max_pos2',
        'confidence_interval_chr2', 'confidence_interval_start2',
        'confidence_interval_end2', 'type', 'score', 'strands']

    order = order + list(set(columns) - set(order))

    df = df[order]

    return df


def write_to_csv(df, filepath):
    df.to_csv(filepath, index=False, na_rep='NaN',
              compression=helpers.get_compression_type_pandas(filepath))


def generate_secondary_table(parsed_data):
    counts = {}

    for brk_call, brk_data in parsed_data:
        breakpoint_id = brk_call[6]

        for dval in brk_data:
            dval = dval.strip().split()

            assert len(dval) == 13

            cell_id = dval[0].split(':')[0]

            if (breakpoint_id, cell_id) not in counts:
                counts[(breakpoint_id, cell_id)] = 0

            counts[(breakpoint_id, cell_id)] += 1

    data = []
    for (brkpt_id, cell_id), count in counts.items():
        data.append((brkpt_id, cell_id, count))

    df = pd.DataFrame(data)

    if df.empty:
        df = pd.DataFrame(columns=['breakpoint_id', 'cell_id', 'count'])
    else:
        df.columns = ['breakpoint_id', 'cell_id', 'count']

    return df


def parse_lumpy(infile, breakpoints, evidence):
    data = group_lumpy_data(infile)
    df = generate_primary_table(data)
    write_to_csv(df, breakpoints)

    data = group_lumpy_data(infile)
    df = generate_secondary_table(data)
    write_to_csv(df, evidence)
