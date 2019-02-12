import pandas as pd


def group_lumpy_data(infile):
    parsed_data = {}

    with open(infile) as data:
        for line in data:
            if not line.startswith("\t"):
                if parsed_data:
                    yield parsed_data
                parsed_data = {}
                key = tuple(line.split()[:-1])
                assert key not in parsed_data
                parsed_data[key] = []
            else:
                parsed_data[key].append(line)


def generate_primary_table(parsed_data):
    data = []

    for lumpy_call in parsed_data:
        assert len(lumpy_call.keys()) == 1

        data.append(lumpy_call.keys()[0])

    df = pd.DataFrame(data)
    df.columns = [
        'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'score',
        'strand1', 'strand2', 'typ', 'ids', 'strands', 'maxv']

    return df


def write_to_h5(df, filepath, tablename):
    with pd.HDFStore(filepath, 'a', complevel=9, complib='blosc') as outfile:
        outfile.put(tablename, df, format='table')


def generate_secondary_table(parsed_data):
    df = []

    for lumpy_call in parsed_data:
        assert len(lumpy_call.keys()) == 1
        key = lumpy_call.keys()[0]
        data = lumpy_call[key]
        refchrom1, refstart1, refend1, refchrom2, refstart2, refend2 = key[:6]


        for dval in data:
            dval = dval.strip().split()

            if len(dval) == 12:
                read_id = "PE"
                chrom1, start1, end1, chrom2, start2, end2 = dval[:6]
            elif len(dval) == 13:
                read_id, chrom1, start1, end1, chrom2, start2, end2 = dval[:7]
            else:
                raise Exception("unknown_format")

            df.append(
                [refchrom1, refstart1, refend1, refchrom2, refstart2, refend2,
                 read_id, chrom1, start1, end1, chrom2, start2, end2]
            )

    df = pd.DataFrame(df)
    df.columns = ['reference_chrom1', 'reference_start1', 'reference_end1', 'reference_chrom2',
                  'reference_start2', 'reference_end2', 'read_id',
                  'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
    return df


def parse_lumpy(infile, outfile):
    data = group_lumpy_data(infile)
    df = generate_primary_table(data)
    write_to_h5(df, outfile, "/lumpy/breakpoints")

    data = group_lumpy_data(infile)
    df = generate_secondary_table(data)
    write_to_h5(df, outfile, "/lumpy/evidence")
