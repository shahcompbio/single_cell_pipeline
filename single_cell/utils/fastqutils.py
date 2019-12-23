import re
from collections import defaultdict
from itertools import islice

from single_cell.utils import helpers


class FastqReader(object):

    def __init__(self, filepath):
        self.file_path = filepath

    def get_read_iterator(self):
        with helpers.getFileHandle(self.file_path) as fq_reader:
            while True:
                fastq_read = list(islice(fq_reader, 4))

                fastq_read = [line for line in fastq_read]

                if not fastq_read:
                    break

                assert len(fastq_read) == 4, 'fastq file format error'

                if not fastq_read[0].startswith('@'):
                    raise ValueError('Expected @ as first character of read name')

                if not fastq_read[2].startswith('+'):
                    raise ValueError('Expected = as first character of read comment')

                yield fastq_read


def _get_read_name(fastq_line1):
    return re.split('/| |\t|#FQST:', fastq_line1.split()[0])[0].rstrip()


class PairedFastqReader(object):
    def __init__(self, r1_path, r2_path):
        self.reader_r1 = FastqReader(r1_path)
        self.reader_r2 = FastqReader(r2_path)

    def get_read_pair_iterator(self):
        read_r1_iter = self.reader_r1.get_read_iterator()
        read_r2_iter = self.reader_r2.get_read_iterator()

        for read_r1, read_r2 in zip(read_r1_iter, read_r2_iter):

            if not read_r1 and not read_r2:
                break

            if not read_r1 or not read_r2:
                raise Exception('mismatching number of reads in R1 and R2')

            assert _get_read_name(read_r1[0]) == _get_read_name(read_r2[0])

            yield read_r1, read_r2


class TaggedFastqReader(FastqReader):
    def __init__(self, fastq_path):
        super(TaggedFastqReader, self).__init__(fastq_path)
        self.indices = None

    def get_read_tag(self, fastq_read):

        read_id = fastq_read[0]
        fq_tag = read_id[read_id.index('FQST:'):]
        fq_tag = fq_tag.strip().split(':')

        if not self.indices:
            if len(fq_tag) > 2:
                self.indices = {i: v for i, v in enumerate(fq_tag[1:-1])}
            else:
                raise Exception('First line in fastq file should have filter explanation')

        flag = map(int, list(fq_tag[-1]))

        flag_map = {self.indices[i]: v for i, v in enumerate(flag)}

        return flag_map

    def add_tag_to_read_comment(self, read, tag=None):
        read_name = _get_read_name(read[0])

        if not tag:
            tag = self.get_read_tag(read)

        tag = ['{}_{}'.format(k, v) for k, v in tag.items()]
        tag = ','.join(tag)
        tag = 'FS:Z:'+tag

        comment = tag

        read[0] = read_name + '\t' + comment + '\n'

        return read

    def filter_read_iterator(self, reference):
        for read in self.get_read_iterator():
            read_tags = self.get_read_tag(read)

            # skip if read maps to multiple genomes
            if len([v for v in read_tags.values() if v]) > 1:
                continue

            if read_tags[reference]:
                yield read
                continue

    def gather_counts(self):
        key_order = None
        counts = defaultdict(int)

        for read in self.get_read_iterator():

            read_tags = self.get_read_tag(read)
            if not key_order:
                key_order = sorted(read_tags.keys())

            flags = zip(key_order, [read_tags[key] for key in key_order])

            counts[flags] += 1

        return counts


class PairedTaggedFastqReader(PairedFastqReader, TaggedFastqReader):
    def __init__(self, fastq_r1, fastq_r2):
        super(PairedTaggedFastqReader, self).__init__(fastq_r1, fastq_r2)
        self.indices = None

    def filter_read_iterator(self, reference):
        for read_1, read_2 in self.get_read_pair_iterator():

            tags_r1 = self.get_read_tag(read_1)
            tags_r2 = self.get_read_tag(read_2)

            # skip if doesnt match
            if not tags_r1[reference] and not tags_r2[reference]:
                continue

            # skip if read maps to multiple genomes
            if len([v for v in tags_r1.values() if v]) > 1:
                continue
            if len([v for v in tags_r2.values() if v]) > 1:
                continue

            r1_nomatch = set(tags_r1.values()) == {0}
            r2_nomatch = set(tags_r2.values()) == {0}

            if tags_r1[reference] and tags_r2[reference]:
                yield read_1, read_2
            elif tags_r1[reference] and r2_nomatch:
                yield read_1, read_2
            elif r1_nomatch and tags_r2[reference]:
                yield read_1, read_2

    def gather_counts(self):
        key_order = None
        counts = {'R1': defaultdict(int), 'R2': defaultdict(int)}

        for read_1, read_2 in self.get_read_pair_iterator():
            tags_r1 = self.get_read_tag(read_1)
            tags_r2 = self.get_read_tag(read_2)

            if not key_order:
                key_order = sorted(tags_r1.keys())

            r1_flags = tuple(zip(key_order, [tags_r1[key] for key in key_order]))
            r2_flags = tuple(zip(key_order, [tags_r2[key] for key in key_order]))

            counts["R1"][r1_flags] += 1
            counts["R2"][r2_flags] += 1

        return counts
