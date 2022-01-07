from collections import defaultdict

import pysam
import yaml
from single_cell.utils import csvutils
from single_cell.workflows.align.dtypes import dtypes


class CoverageMetrics(object):
    def __init__(
            self,
            bamfile,
            filter_unpaired=False,
            filter_duplicates=False,
            filter_supplementary=False,
            filter_secondary=False,
            min_mapping_qual=0,
            min_base_qual=0
    ):
        self.bamfile = bamfile
        self.filter_unpaired = filter_unpaired
        self.filter_duplicates = filter_duplicates
        self.filter_supplementary = filter_supplementary
        self.filter_secondary = filter_secondary
        self.min_mapping_qual = min_mapping_qual
        self.min_base_qual = min_base_qual
        self._bam_reader = None

    def __enter__(self):
        self._bam_reader = self._get_bam_reader()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._bam_reader.close()

    def _get_bam_reader(self):
        return pysam.AlignmentFile(self.bamfile, 'rb')

    @property
    def genome_length(self):
        lengths = [val['LN'] for val in self._bam_reader.header['SQ']]
        return sum(lengths)

    def _get_read_intervals(self, read):
        if self.min_base_qual == 0:
            return [(read.reference_start, read.reference_end)]

        regions = []
        start = None

        read_quals = read.query_alignment_qualities
        if read.is_reverse:
            read_quals = read_quals[::-1]
        for i, qual in zip(range(read.reference_start, read.reference_end), read_quals):
            if start is None and qual >= self.min_base_qual:
                start = i
            else:
                if start is not None and qual < self.min_base_qual:
                    regions.append((start, i))
                    start = None
        if start is not None:
            regions.append((start, i))

        assert len(regions) > 0
        return regions

    def _filter_reads(self, read):
        if not read.is_paired and self.filter_unpaired:
            return True

        if read.is_duplicate and self.filter_duplicates:
            return True

        if read.is_supplementary and self.filter_supplementary:
            return True

        if read.is_secondary and self.filter_secondary:
            return True

        if read.mapping_quality < self.min_mapping_qual:
            return True

        # read is completely unmapped
        if read.reference_end is None:
            return True

        return False

    @staticmethod
    def _merge_overlapping_intervals(coords):
        coords.sort(key=lambda x: x[0])

        if len(coords) < 2:
            return coords
        elif len(coords) == 2:
            if coords[1][0] <= coords[0][1]:
                return [[coords[0][0], max(coords[0][1], coords[1][1])]]
            else:
                return coords
        else:
            i = 0
            while i < len(coords) - 1:
                if i < len(coords) - 1:
                    if coords[i + 1][0] <= coords[i][1]:
                        start = coords[i][0]
                        end = max(coords[i][1], coords[i + 1][1])
                        del coords[i + 1]
                        coords[i] = [start, end]
                    else:
                        i += 1
                else:
                    i += 1
            return coords

    def generate_data(self):
        read_dict = defaultdict(lambda: defaultdict(list))
        for read in self._bam_reader.fetch():
            if self._filter_reads(read) is True:
                continue

            regions = self._get_read_intervals(read)
            read_dict[read.query_name][read.reference_name].extend(regions)
        return read_dict

    def get_coverage(self, read_dict):
        total_length = 0
        for rname, chromdata in read_dict.items():
            for chrom, regions in chromdata.items():
                regions = self._merge_overlapping_intervals(regions)
                chrom_length = sum([v[1] - v[0] for v in regions])

                total_length += chrom_length

        return float(total_length) / self.genome_length

    def main(self):
        if self._bam_reader is None:
            self._bam_reader = self._get_bam_reader()

        read_data = self.generate_data()
        coverage = self.get_coverage(read_data)
        return coverage


def expected_and_aligned_coverage(bamfile):
    with CoverageMetrics(bamfile) as cov:
        genome_length = cov.genome_length

    bam = pysam.AlignmentFile(bamfile, 'rb')
    expected_length = 0
    aligned_length = 0
    for read in bam:
        expected_length += read.query_length
        aligned_length += 0 if read.reference_length is None else read.reference_length

    return expected_length / genome_length, aligned_length / genome_length


def get_coverage_data(bamfile, output, cell_id, mapping_qual=10, base_qual=10):
    outdata = {'cell_id': cell_id}

    expected_cov, aligned_cov = expected_and_aligned_coverage(bamfile)

    outdata['expected'] = expected_cov
    outdata['aligned'] = aligned_cov

    with CoverageMetrics(bamfile) as cov:
        outdata['overlap_with_dups'] = cov.main()

    with CoverageMetrics(bamfile, filter_duplicates=True) as cov:
        outdata['overlap_without_dups'] = cov.main()

    with CoverageMetrics(
            bamfile, filter_duplicates=True, filter_secondary=True,
            filter_supplementary=True, filter_unpaired=True
    ) as cov:
        outdata['overlap_with_all_filters'] = cov.main()

    with CoverageMetrics(
            bamfile, filter_duplicates=True, filter_secondary=True,
            filter_supplementary=True, filter_unpaired=True,
            min_base_qual=base_qual, min_mapping_qual=mapping_qual
    ) as cov:
        outdata['overlap_with_all_filters_and_qual'] = cov.main()

    with open(output, 'wt') as writer:
        yaml.dump(outdata, writer)


def annotate_coverage_metrics(metrics, coverage_yaml, output):
    data = {}

    for cell_id, filename in coverage_yaml.items():
        with open(filename, 'rt') as reader:
            covdata = yaml.load(reader)
            if 'cell_id' in covdata:
                assert covdata['cell_id'] == cell_id
                del covdata['cell_id']
            data[cell_id] = covdata

    csvutils.annotate_csv(
        metrics, data, output, dtypes()['metrics']
    )
