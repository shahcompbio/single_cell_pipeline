import gzip
import os

from single_cell.utils import fastqutils
from single_cell.utils import helpers


def get_basename(filepath):
    filepath_base = os.path.basename(filepath)

    if filepath_base.endswith('.fastq.gz'):
        filepath_base = filepath_base[:-len('.fastq.gz')]
        extension = '.fastq.gz'
    elif filepath_base.endswith('.fq.gz'):
        filepath_base = filepath_base[:-len('.fq.gz')]
        extension = '.fastq.gz'
    elif filepath_base.endswith('.fastq'):
        filepath_base = filepath_base[:-len('.fastq')]
        extension = '.fastq'
    elif filepath_base.endswith('.fq'):
        filepath_base = filepath_base[:-len('.fq')]
        extension = '.fastq'
    else:
        raise Exception('unknown file format. {}'.format(filepath))
    return filepath_base, extension


def regroup_needed(params):
    for genome in params['genomes']:
        genome_paths = genome['paths']
        if isinstance(genome_paths, list) and len(genome_paths) > 1:
            return True
    return False


def generate_fastqscreen_config(filepath, params):
    with open(filepath, 'w') as config_writer:
        for genome in params['genomes']:
            genome_name = genome['name']
            if '_' in genome_name:
                raise Exception('_ not allowed in genome name in fastqscreen config')

            genome_paths = genome['paths']

            if isinstance(genome_paths, list) and len(genome_paths) > 1:
                for i, genome_path in enumerate(genome_paths):
                    outstr = '\t'.join(['DATABASE', f'{genome_name}_path{i}', genome_path]) + '\n'
                    config_writer.write(outstr)
            else:
                config_writer.write('DATABASE\t{}\t{}\n'.format(genome_name, genome_paths))


def update_read_tags(read, newtags):
    start, tag = read[0].split('#FQST:')

    if len(tag.split(':')) > 1:
        tag = ':'.join(newtags.keys()) + ':' + ''.join([str(v) for v in newtags.values()])
    else:
        tag = ''.join([str(v) for v in newtags.values()])

    read[0] = '{}#FQST:{}\n'.format(start, tag)

    return read


def regroup_genomes(input_fastq, output_fastq):
    reader = fastqutils.TaggedFastqReader(input_fastq)

    assert output_fastq.endswith('gz')
    output = gzip.open(output_fastq, 'wt')

    for read in reader.get_read_iterator():
        tags = reader.get_read_tag(read)

        newtags = {}
        for tag, val in tags.items():

            if '_path' in tag:
                tag = tag.split('_path')[0]

            if tag in newtags:
                newtags[tag] = max(val, newtags[tag])
            else:
                newtags[tag] = val

        read = update_read_tags(read, newtags)

        for line in read:
            output.write(line)


def filter_tag_reads(
        input_r1, input_r2, output_r1, output_r2, params
):
    genomes = [v['name'] for v in params['genomes']]

    if not params['filter_tags']:
        filter_tags = set()
    else:
        filter_tags = set(params['filter_tags'])

    reader = fastqutils.PairedTaggedFastqReader(input_r1, input_r2)

    with helpers.getFileHandle(output_r1, 'wt') as writer_r1, helpers.getFileHandle(output_r2, 'wt') as writer_r2:
        for read_1, read_2 in reader.filter_read_iterator(genomes, filter_tags):

            read_1 = reader.add_tag_to_read_comment(read_1)
            read_2 = reader.add_tag_to_read_comment(read_2)

            for line in read_1:
                writer_r1.write(line)

            for line in read_2:
                writer_r2.write(line)
