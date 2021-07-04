import os
import shutil
from collections import defaultdict

import pandas as pd
import pypeliner
import single_cell.workflows.align.fastqscreen_utils as utils
from single_cell.utils import csvutils
from single_cell.utils import fastqutils
from single_cell.utils import helpers
from single_cell.workflows.align.dtypes import fastqscreen_dtypes


def merge_fastq_screen_counts(
        all_detailed_counts, all_summary_counts,
        merged_detailed_counts, merged_summary_counts,
        fastqscreen_config
):
    genome_labels = [genome['name'] for genome in fastqscreen_config['genomes']]

    all_detailed_counts = helpers.flatten(all_detailed_counts)
    all_detailed_counts = [pd.read_csv(file) for file in all_detailed_counts if not helpers.is_empty(file)]
    df = pd.concat(all_detailed_counts)

    index_cols = [v for v in df.columns.values if v != "count"]

    df['count'] = df.groupby(index_cols)['count'].transform('sum')

    df = df.drop_duplicates(subset=index_cols)

    csvutils.write_dataframe_to_csv_and_yaml(
        df, merged_detailed_counts,
        fastqscreen_dtypes(genome_labels)['fastqscreen_detailed'], write_header=True
    )

    all_summary_counts = helpers.flatten(all_summary_counts)
    all_summary_counts = [pd.read_csv(file) for file in all_summary_counts if not helpers.is_empty(file)]
    df = pd.concat(all_summary_counts)

    update_cols = [v for v in df.columns.values if v != 'cell_id']

    for colname in update_cols:
        df[colname] = df.groupby('cell_id')[colname].transform('sum')

    df = df.drop_duplicates(subset=['cell_id'])

    csvutils.write_dataframe_to_csv_and_yaml(
        df, merged_summary_counts,
        fastqscreen_dtypes(genome_labels)['metrics'], write_header=True
    )


def run_fastq_screen_paired_end(fastq_r1, fastq_r2, tempdir, params):
    r1_basename, r1_ext = utils.get_basename(fastq_r1)
    tagged_fastq_r1 = os.path.join(tempdir, '{}.tagged{}'.format(r1_basename, r1_ext))

    r2_basename, r2_ext = utils.get_basename(fastq_r2)
    tagged_fastq_r2 = os.path.join(tempdir, '{}.tagged{}'.format(r2_basename, r2_ext))

    if helpers.is_empty(fastq_r1):
        shutil.copy(fastq_r1, tagged_fastq_r1)
        shutil.copy(fastq_r2, tagged_fastq_r2)
        return tagged_fastq_r1, tagged_fastq_r2

    config = os.path.join(tempdir, 'fastq_screen.config')
    utils.generate_fastqscreen_config(config, params)

    pypeliner.commandline.execute(
        'fastq_screen',
        '--aligner', params['aligner'],
        '--conf', config,
        '--outdir', tempdir,
        '--tag',
        fastq_r1, fastq_r2,
    )

    if utils.regroup_needed(params):
        fixed_fastq_r1 = os.path.join(tempdir, '{}.tagged.fixed{}'.format(r1_basename, r1_ext))
        fixed_fastq_r2 = os.path.join(tempdir, '{}.tagged.fixed{}'.format(r2_basename, r2_ext))
        utils.regroup_genomes(tagged_fastq_r1, fixed_fastq_r1)
        utils.regroup_genomes(tagged_fastq_r2, fixed_fastq_r2)
        return fixed_fastq_r1, fixed_fastq_r2
    else:
        return tagged_fastq_r1, tagged_fastq_r2


def write_detailed_counts(counts, outfile, cell_id, fastqscreen_params):
    header = None

    genomes = [genome['name'] for genome in fastqscreen_params['genomes']]

    with helpers.getFileHandle(outfile, 'wt') as writer:

        for read_end, read_end_counts in counts.items():

            if not read_end_counts and not header:
                outstr = ['cell_id', 'readend'] + genomes + ['count']
                writer.write(','.join(outstr) + '\n')
                header = 1
                continue

            if not header:
                outstr = ['cell_id', 'readend']
                outstr += [v[0] for v in list(read_end_counts.keys())[0]]
                outstr += ['count']
                writer.write(','.join(outstr) + '\n')
                header = 1

            for flags, count in read_end_counts.items():
                outstr = [cell_id, read_end]
                outstr += [v[1] for v in flags]
                outstr += [count]
                writer.write(','.join(map(str, outstr)) + '\n')


def write_summary_counts(counts, outfile, cell_id, fastqscreen_params):
    genomes = [genome['name'] for genome in fastqscreen_params['genomes']]

    summary_counts = defaultdict(int)
    for read_end, read_end_counts in counts.items():
        for flags, count in read_end_counts.items():
            hit_orgs = [v[0] for v in flags if v[1] > 0]

            for org in hit_orgs:
                summary_counts[org] += count

            if len(hit_orgs) > 1:
                for org in hit_orgs:
                    summary_counts['{}_multihit'.format(org)] += count
            elif len(hit_orgs) == 0:
                summary_counts['nohit'] += count

    with helpers.getFileHandle(outfile, 'wt') as writer:
        if not summary_counts:
            columns = ['cell_id']
            columns += ['fastqscreen_' + genome for genome in genomes]
            columns += ['fastqscreen_nohit']
            header = ','.join(columns) + '\n'
            writer.write(header)
            data = [0] * len(columns)
            data[0] = cell_id
            data = [str(v) for v in data]
            data = ','.join(data) + '\n'
            writer.write(data)
            return

        keys = sorted(summary_counts.keys())
        header = ['cell_id'] + ['fastqscreen_{}'.format(key) for key in keys]
        header = ','.join(header) + '\n'
        writer.write(header)

        values = [cell_id] + [summary_counts[v] for v in keys]
        values = ','.join(map(str, values)) + '\n'
        writer.write(values)




def organism_filter(
        fastq_r1, fastq_r2, filtered_fastq_r1, filtered_fastq_r2,
        detailed_metrics, summary_metrics, tempdir, cell_id, params
):
    # fastq screen tries to skip if files from old runs are available
    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)

    helpers.makedirs(tempdir)

    tagged_fastq_r1, tagged_fastq_r2 = run_fastq_screen_paired_end(
        fastq_r1, fastq_r2, tempdir, params,
    )

    reader = fastqutils.PairedTaggedFastqReader(tagged_fastq_r1, tagged_fastq_r2)
    counts = reader.gather_counts()

    write_detailed_counts(counts, detailed_metrics, cell_id, params)
    write_summary_counts(counts, summary_metrics, cell_id, params)

    utils.filter_tag_reads(
        tagged_fastq_r1, tagged_fastq_r2, filtered_fastq_r1,
        filtered_fastq_r2, params
    )
