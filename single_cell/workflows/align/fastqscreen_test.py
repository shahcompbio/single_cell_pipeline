import gzip
import os

from single_cell.workflows.align.fastqscreen import filter_tag_reads


def simulate_paired_fastq(r1, r2, flags, values, id):
    flags = ':'.join(flags)
    vals = ''.join(map(str, values))

    readid = '@HISEQ{}_144:5:1101:6674:43220#FQST:{}:{}'.format(id, flags, vals)

    seq1 = 'ACGT' * 30

    seq2 = 'GCTT' * 30

    r1.write(readid + '\n')
    r1.write(seq1 + '\n')
    r1.write('+\n')
    r1.write(seq2 + '\n')

    r2.write(readid + '\n')
    r2.write(seq1 + '\n')
    r2.write('+\n')
    r2.write(seq2 + '\n')


def get_line_count(filepath):
    count = 0
    with gzip.open(filepath, 'rt') as reader:
        for _ in reader:
            count += 1
    return count


def test_no_filter(tempdir):
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    r1 = os.path.join(tempdir, 'R1.fastq.gz')
    r2 = os.path.join(tempdir, 'R2.fastq.gz')

    r1_filt = os.path.join(tempdir, 'R1_filtered.fastq.gz')
    r2_filt = os.path.join(tempdir, 'R2_filtered.fastq.gz')

    with gzip.open(r1, 'wt') as r1_out, gzip.open(r2, 'wt') as r2_out:
        simulate_paired_fastq(r1_out, r2_out, ['grch37', 'mm10', 'salmon', 'bact'], [1, 0, 0, 0], '101')

    filter_tag_reads(
        r1, r2, r1_filt, r2_filt,
        {'grch37': False, 'mm10': False, 'salmon': False, 'bact': False},
        {'grch37': False, 'mm10': False, 'salmon': False, 'bact': False}
    )

    assert get_line_count(r1_filt) == get_line_count(r1)
    assert get_line_count(r2_filt) == get_line_count(r2)


def test_filter_out_bact_inclusive(tempdir):
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    r1 = os.path.join(tempdir, 'R1.fastq.gz')
    r2 = os.path.join(tempdir, 'R2.fastq.gz')

    r1_filt = os.path.join(tempdir, 'R1_filtered.fastq.gz')
    r2_filt = os.path.join(tempdir, 'R2_filtered.fastq.gz')

    with gzip.open(r1, 'wt') as r1_out, gzip.open(r2, 'wt') as r2_out:
        simulate_paired_fastq(r1_out, r2_out, ['grch37', 'mm10', 'salmon', 'bact'], [0, 0, 0, 1], '101')
        simulate_paired_fastq(r1_out, r2_out, ['grch37', 'mm10', 'salmon', 'bact'], [1, 0, 0, 0], '103')
        simulate_paired_fastq(r1_out, r2_out, ['grch37', 'mm10', 'salmon', 'bact'], [1, 0, 0, 1], '104')

    filter_tag_reads(
        r1, r2, r1_filt, r2_filt,
        {'grch37': False, 'mm10': False, 'salmon': False, 'bact': True},
        {'grch37': False, 'mm10': False, 'salmon': False, 'bact': False}
    )

    assert get_line_count(r1_filt) == 4
    assert get_line_count(r2_filt) == 4


def test_filter_out_bact_exclusive(tempdir):
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    r1 = os.path.join(tempdir, 'R1.fastq.gz')
    r2 = os.path.join(tempdir, 'R2.fastq.gz')

    r1_filt = os.path.join(tempdir, 'R1_filtered.fastq.gz')
    r2_filt = os.path.join(tempdir, 'R2_filtered.fastq.gz')

    with gzip.open(r1, 'wt') as r1_out, gzip.open(r2, 'wt') as r2_out:
        simulate_paired_fastq(r1_out, r2_out, ['grch37', 'mm10', 'salmon', 'bact'], [0, 0, 0, 1], '101')
        simulate_paired_fastq(r1_out, r2_out, ['grch37', 'mm10', 'salmon', 'bact'], [1, 0, 0, 0], '103')
        simulate_paired_fastq(r1_out, r2_out, ['grch37', 'mm10', 'salmon', 'bact'], [1, 0, 0, 1], '104')

    filter_tag_reads(
        r1, r2, r1_filt, r2_filt,
        {'grch37': False, 'mm10': False, 'salmon': False, 'bact': False},
        {'grch37': False, 'mm10': False, 'salmon': False, 'bact': True}
    )

    assert get_line_count(r1_filt) == 8
    assert get_line_count(r2_filt) == 8


def test_filter_out_bact_inclusive_multi(tempdir):
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    r1 = os.path.join(tempdir, 'R1.fastq.gz')
    r2 = os.path.join(tempdir, 'R2.fastq.gz')

    r1_filt = os.path.join(tempdir, 'R1_filtered.fastq.gz')
    r2_filt = os.path.join(tempdir, 'R2_filtered.fastq.gz')

    with gzip.open(r1, 'wt') as r1_out, gzip.open(r2, 'wt') as r2_out:
        simulate_paired_fastq(r1_out, r2_out, ['grch37', 'mm10', 'salmon', 'bact'], [0, 0, 0, 1], '101')
        simulate_paired_fastq(r1_out, r2_out, ['grch37', 'mm10', 'salmon', 'bact'], [0, 0, 1, 0], '102')
        simulate_paired_fastq(r1_out, r2_out, ['grch37', 'mm10', 'salmon', 'bact'], [1, 0, 0, 0], '103')
        simulate_paired_fastq(r1_out, r2_out, ['grch37', 'mm10', 'salmon', 'bact'], [1, 0, 0, 1], '104')
        simulate_paired_fastq(r1_out, r2_out, ['grch37', 'mm10', 'salmon', 'bact'], [1, 0, 1, 1], '105')

    filter_tag_reads(
        r1, r2, r1_filt, r2_filt,
        {'grch37': False, 'mm10': False, 'salmon': True, 'bact': True},
        {'grch37': False, 'mm10': False, 'salmon': False, 'bact': False}
    )

    assert get_line_count(r1_filt) == 4
    assert get_line_count(r2_filt) == 4


def test_filter_out_bact_exclusive_multi(tempdir):
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    r1 = os.path.join(tempdir, 'R1.fastq.gz')
    r2 = os.path.join(tempdir, 'R2.fastq.gz')

    r1_filt = os.path.join(tempdir, 'R1_filtered.fastq.gz')
    r2_filt = os.path.join(tempdir, 'R2_filtered.fastq.gz')

    with gzip.open(r1, 'wt') as r1_out, gzip.open(r2, 'wt') as r2_out:
        simulate_paired_fastq(r1_out, r2_out, ['grch37', 'mm10', 'salmon', 'bact'], [0, 0, 0, 1], '101')
        simulate_paired_fastq(r1_out, r2_out, ['grch37', 'mm10', 'salmon', 'bact'], [0, 0, 1, 0], '102')
        simulate_paired_fastq(r1_out, r2_out, ['grch37', 'mm10', 'salmon', 'bact'], [1, 0, 0, 1], '103')
        simulate_paired_fastq(r1_out, r2_out, ['grch37', 'mm10', 'salmon', 'bact'], [1, 0, 1, 0], '104')
        simulate_paired_fastq(r1_out, r2_out, ['grch37', 'mm10', 'salmon', 'bact'], [0, 0, 0, 0], '105')

    filter_tag_reads(
        r1, r2, r1_filt, r2_filt,
        {'grch37': False, 'mm10': False, 'salmon': False, 'bact': False},
        {'grch37': False, 'mm10': False, 'salmon': True, 'bact': True}
    )

    assert get_line_count(r1_filt) == 12
    assert get_line_count(r2_filt) == 12


tempdir = 'temp'

test_no_filter(os.path.join(tempdir, 'nofilter'))
test_filter_out_bact_inclusive(os.path.join(tempdir, 'bact_inclusive'))
test_filter_out_bact_exclusive(os.path.join(tempdir, 'bact_exclusive'))
test_filter_out_bact_inclusive_multi(os.path.join(tempdir, 'bact_exclusive_multu'))
test_filter_out_bact_exclusive_multi(os.path.join(tempdir, 'bact_exclusive_multi'))
