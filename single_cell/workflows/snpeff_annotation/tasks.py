import os
import re
from collections import OrderedDict

import pandas as pd
import pypeliner
import vcf
from single_cell.utils import csvutils
from single_cell.workflows.snpeff_annotation.dtypes import dtypes


def run_snpeff(db, data_dir, in_vcf_file, out_file, classic_mode=True):
    os.environ['MALLOC_ARENA_MAX'] = '2'
    data_dir = os.path.abspath(data_dir)

    cmd = [
        'snpEff',
        '-noStats',
        '-noLog',
        '-Xms2g',
        '-Xmx5g',
        '-hgvs1LetterAa',
        '-dataDir',
        data_dir,
    ]

    if classic_mode:
        cmd.append('-classic')

    cmd.extend([
        db,
        in_vcf_file,
        '>',
        out_file
    ])

    pypeliner.commandline.execute(*cmd)


class ClassicSnpEffParser(object):

    def __init__(self, file_name):
        self._reader = vcf.Reader(filename=file_name)

        self.fields = self._get_field_names()

        self._buffer = []

        self._effect_matcher = re.compile(r'(.*)\(')

        self._fields_matcher = re.compile(r'\((.*)\)')

    def __iter__(self):
        while True:
            try:
                yield self.next()
            except StopIteration:
                break

    def next(self):
        while len(self._buffer) == 0:
            record = next(self._reader)

            if 'EFF' not in record.INFO:
                continue

            for row in self._parse_record(record):
                self._buffer.append(row)

        return self._buffer.pop(0)

    def _get_field_names(self):
        fields = []

        match = re.search(r'\((.*)\[', self._reader.infos['EFF'].desc)

        for x in match.groups()[0].split('|'):
            fields.append(x.strip().lower())

        return fields

    def _parse_record(self, record):
        for annotation in record.INFO['EFF']:
            effect = self._effect_matcher.search(annotation).groups()[0]

            out_row = OrderedDict((
                ('chrom', record.CHROM),
                ('coord', record.POS),
                ('ref', record.REF),
                ('alt', ','.join([str(x) for x in record.ALT])),
                ('effect', effect),
            ))

            fields = self._fields_matcher.search(annotation).groups()[0].split('|')

            for i, key in enumerate(self.fields):
                out_row[key] = fields[i]

            yield out_row


def convert_vcf_to_table(in_file, out_file):
    data = []

    parser = ClassicSnpEffParser(in_file)

    for row in parser:
        data.append(row)

    data = pd.DataFrame(data)

    csvutils.write_dataframe_to_csv_and_yaml(data, out_file, dtypes())
