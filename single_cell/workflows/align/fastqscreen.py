import os
import shutil
from collections import defaultdict

import pandas as pd
import pypeliner
from single_cell.utils import csvutils
from single_cell.utils import fastqutils
from single_cell.utils import helpers
from single_cell.workflows.align.dtypes import dtypes


def merge_fastq_screen_counts(
        all_detailed_counts, all_summary_counts, merged_detailed_counts, merged_summary_counts
):
    if isinstance(all_detailed_counts, dict):
        all_detailed_counts = all_detailed_counts.values()

    detailed_data = []
    for countsfile in all_detailed_counts:
        if os.stat(countsfile).st_size == 0:
            continue
        detailed_data.append(pd.read_csv(countsfile))

    df = pd.concat(detailed_data)

    index_cols = [v for v in df.columns.values if v != "count"]

    df['count'] = df.groupby(index_cols)['count'].transform('sum')

    df = df.drop_duplicates(subset=index_cols)

    csvutils.write_dataframe_to_csv_and_yaml(
        df, merged_detailed_counts, dtypes()['fastqscreen_detailed'], write_header=True
    )

    if isinstance(all_summary_counts, dict):
        all_summary_counts = all_summary_counts.values()

    summary_counts = [pd.read_csv(countsfile) for countsfile in all_summary_counts]

    df = pd.concat(summary_counts)

    update_cols = [v for v in df.columns.values if v != 'cell_id']

    for colname in update_cols:
        df[colname] = df.groupby('cell_id')[colname].transform('sum')

    df = df.drop_duplicates(subset=['cell_id'])

    csvutils.write_dataframe_to_csv_and_yaml(
        df, merged_summary_counts, dtypes()['metrics'], write_header=True
    )


def run_fastq_screen_paired_end(fastq_r1, fastq_r2, tempdir, params, docker_image=None):
    def get_basename(filepath):
        filepath_base = os.path.basename(filepath)

        if filepath_base.endswith('.fastq.gz'):
            filepath_base = filepath_base[:-len('.fastq.gz')]
        elif filepath_base.endswith('.fq.gz'):
            filepath_base = filepath_base[:-len('.fq.gz')]
        elif filepath_base.endswith('.fastq'):
            filepath_base = filepath_base[:-len('.fastq')]
        elif filepath_base.endswith('.fq'):
            filepath_base = filepath_base[:-len('.fq')]
        else:
            raise Exception('unknown file format. {}'.format(filepath))
        return filepath_base

    basename = get_basename(fastq_r1)
    tagged_fastq_r1 = os.path.join(tempdir, '{}.tagged.fastq.gz'.format(basename))

    basename = get_basename(fastq_r2)
    tagged_fastq_r2 = os.path.join(tempdir, '{}.tagged.fastq.gz'.format(basename))

    # fastq screen fails if run on empty files
    with helpers.getFileHandle(fastq_r1) as reader:
        if not reader.readline():
            shutil.copy(fastq_r1, tagged_fastq_r1)
            shutil.copy(fastq_r2, tagged_fastq_r2)
            return tagged_fastq_r1, tagged_fastq_r2

    config = os.path.join(tempdir, 'fastq_screen.config')

    with open(config, 'w') as config_writer:
        for genome in params['genomes']:
            genome_name = genome['name']
            genome_path = genome['path']
            outstr = '\t'.join(['DATABASE', genome_name, genome_path]) + '\n'
            config_writer.write(outstr)

    cmd = [
        'fastq_screen',
        '--aligner', params['aligner'],
        '--conf', config,
        '--outdir', tempdir,
        '--tag',
        fastq_r1,
        fastq_r2,
    ]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    return tagged_fastq_r1, tagged_fastq_r2


def write_detailed_counts(counts, outfile, cell_id, fastqscreen_params):
    header = None

    genomes = [genome['name'] for genome in fastqscreen_params['genomes']]

    with helpers.getFileHandle(outfile, 'wt') as writer:

        for read_end, read_end_counts in counts.items():

            if not read_end_counts:
                outstr = ['cell_id', 'readend'] + genomes + ['count']
                writer.write(','.join(outstr) + '\n')
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
            return

        keys = sorted(summary_counts.keys())
        header = ['cell_id'] + ['fastqscreen_{}'.format(key) for key in keys]
        header = ','.join(header) + '\n'
        writer.write(header)

        values = [cell_id] + [summary_counts[v] for v in keys]
        values = ','.join(map(str, values)) + '\n'
        writer.write(values)


def filter_reads(
        input_r1, input_r2, output_r1, output_r2, reference
):
    reader = fastqutils.PairedTaggedFastqReader(input_r1, input_r2)

    with helpers.getFileHandle(output_r1, 'wt') as writer_r1, helpers.getFileHandle(output_r2, 'wt') as writer_r2:
        for read_1, read_2 in reader.filter_read_iterator(reference):

            read_1 = reader.add_tag_to_read_comment(read_1)
            read_2 = reader.add_tag_to_read_comment(read_2)

            for line in read_1:
                writer_r1.write(line)

            for line in read_2:
                writer_r2.write(line)


def re_tag_reads(infile, outfile):
    reader = fastqutils.TaggedFastqReader(infile)

    with helpers.getFileHandle(outfile, 'wt') as writer:

        for read in reader.get_read_iterator():
            read = reader.add_tag_to_read_comment(read)

            for line in read:
                writer.write(line)


def organism_filter(
        fastq_r1, fastq_r2, filtered_fastq_r1, filtered_fastq_r2,
        detailed_metrics, summary_metrics, tempdir, cell_id, params,
        reference, docker_image=None, filter_contaminated_reads=False,
):
    # fastq screen tries to skip if files from old runs are available
    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)

    helpers.makedirs(tempdir)

    tagged_fastq_r1, tagged_fastq_r2 = run_fastq_screen_paired_end(
        fastq_r1, fastq_r2, tempdir, params, docker_image=docker_image
    )

    reader = fastqutils.PairedTaggedFastqReader(tagged_fastq_r1, tagged_fastq_r2)
    counts = reader.gather_counts()

    write_detailed_counts(counts, detailed_metrics, cell_id, params)
    write_summary_counts(counts, summary_metrics, cell_id, params)

    if filter_contaminated_reads:
        filter_reads(
            tagged_fastq_r1, tagged_fastq_r2, filtered_fastq_r1,
            filtered_fastq_r2, reference
        )
    else:
        # use the full tagged fastq downstream
        # with organism type information in readname
        re_tag_reads(tagged_fastq_r1, filtered_fastq_r1)
        re_tag_reads(tagged_fastq_r2, filtered_fastq_r2)
