'''
Extract metrics table.
'''

from __future__ import division

import os

import pandas as pd

from single_cell.utils import csvutils

class CollectMetrics(object):
    def __init__(
            self, wgs_metrics, insert_metrics, flagstat_metrics,
            markdups_metrics, output, sample_id, dtypes
    ):
        self.wgs_metrics = wgs_metrics
        self.flagstat_metrics = flagstat_metrics
        self.insert_metrics = insert_metrics
        self.markdups_metrics = markdups_metrics
        self.output = output
        self.sample_id = sample_id
        self.dtypes = dtypes

    def extract_wgs_metrics(self):
        """
        get the coverage_depth (mean_coverage column)
        get the coverage_breadth (count/genome_territory)
        """

        mfile = open(self.wgs_metrics)

        metrics = []
        hist = {}

        addmetrics = False
        addhist = False

        for line in mfile:
            if line.strip() == '':
                continue
            if line.startswith('## METRICS CLASS'):
                addmetrics = True
                addhist = False
                continue

            if line.startswith('## HISTOGRAM'):
                addhist = True
                addmetrics = False
                continue

            if addmetrics:
                metrics.append(line.strip().split('\t'))
            if addhist:
                line = line.strip().split('\t')
                if line[0] == 'coverage':
                    continue
                hist[int(line[0])] = int(line[1])

        mfile.close()
        header, data = metrics

        header = [v.lower() for v in header]
        header = {v: i for i, v in enumerate(header)}

        gen_territory = int(data[header['genome_territory']])
        cov_depth = float(data[header['mean_coverage']])
        count = int(hist[0])
        cov_breadth = (gen_territory - count) / gen_territory

        return cov_breadth, cov_depth

    def extract_flagstat_metrics(self):
        """
        extract from flagstat
        """

        df = pd.read_csv(
            self.flagstat_metrics,
            sep=r'\s\+\s0\s',
            header=None,
            names=['value', 'type'],
            engine='python'
        )

        tot_reads = df[df['type'] == 'in total (QC-passed reads + QC-failed reads)']['value']
        tot_mpd_reads = df[
            (df['type'].str.contains('mapped') == True) & (df['type'].str.contains('mate mapped') == False)]
        tot_dup_reads = df[df['type'] == 'duplicates']['value']
        tot_prop_paired = df[df['type'].str.contains('properly paired')]

        assert len(tot_reads) == 1
        assert len(tot_mpd_reads) == 1
        assert len(tot_dup_reads) == 1
        assert len(tot_prop_paired) == 1

        tot_reads = tot_reads.iloc[0]
        tot_mpd_reads = tot_mpd_reads['value'].iloc[0]
        tot_dup_reads = tot_dup_reads.iloc[0]
        tot_prop_paired = tot_prop_paired['value'].iloc[0]

        return tot_reads, tot_mpd_reads, tot_dup_reads, tot_prop_paired

    def extract_duplication_metrics(self):
        """
        extract from markdups
        """

        mfile = open(self.markdups_metrics)

        targetlines = []

        line = mfile.readline()

        while line != '':
            if line.startswith('## METRICS CLASS'):
                targetlines.append(mfile.readline().strip('\n').split('\t'))
                targetlines.append(mfile.readline().strip('\n').split('\t'))
                break
            line = mfile.readline()

        mfile.close()

        header, data = targetlines

        header = [v.lower() for v in header]
        header = {v: i for i, v in enumerate(header)}

        unprd_mpd_rds = int(data[header['unpaired_reads_examined']])
        prd_mpd_rds = int(data[header['read_pairs_examined']])
        unprd_dup_rds = int(data[header['unpaired_read_duplicates']])
        prd_dup_rds = int(data[header['read_pair_duplicates']])
        unmpd_rds = data[header['unmapped_reads']]
        est_lib_size = data[header['estimated_library_size']]

        rd_pair_opt_dup = int(data[header['read_pair_optical_duplicates']])

        try:
            perc_dup_reads = (unprd_dup_rds + ((prd_dup_rds + rd_pair_opt_dup) * 2)) / (
                        unprd_mpd_rds + (prd_mpd_rds * 2))
        except ZeroDivisionError:
            perc_dup_reads = 0

        outdata = (unprd_mpd_rds, prd_mpd_rds, unprd_dup_rds, prd_dup_rds,
                   unmpd_rds, perc_dup_reads, est_lib_size)

        outdata = tuple([0 if val == '' else val for val in outdata])
        return outdata

    def extract_insert_metrics(self):
        ''' Extract median and mean insert size '''

        # picardtools insertmetrics completes with code 0 and doesn't generate metrics file
        # if inputs don't have sufficient read count
        if not os.path.isfile(self.insert_metrics):
            return 0, 0, 0

        # if the insert metrics fails due to low coverage
        if open(self.insert_metrics).readline().startswith("## FAILED"):
            return 0, 0, 0

        mfile = open(self.insert_metrics)

        targetlines = []

        line = mfile.readline()

        while line != '':
            if line.startswith('## METRICS CLASS'):
                targetlines.append(mfile.readline().strip().split('\t'))
                targetlines.append(mfile.readline().strip().split('\t'))
                break
            line = mfile.readline()

        mfile.close()

        header, data = targetlines

        header = [v.lower() for v in header]
        header = {v: i for i, v in enumerate(header)}

        median_ins_size = data[header['median_insert_size']]
        mean_ins_size = data[header['mean_insert_size']]
        std_dev_ins_size = data[header['standard_deviation']]

        median_ins_size = 0 if median_ins_size == '?' else median_ins_size
        mean_ins_size = 0 if mean_ins_size == '?' else mean_ins_size
        std_dev_ins_size = 0 if std_dev_ins_size == '?' else std_dev_ins_size

        return median_ins_size, mean_ins_size, std_dev_ins_size

    def write_data(self, header, data):
        """
        write to the output
        """
        assert len(header) == len(data)
        # replace empty vals with NA
        data = [v if v != '' else 'NA' for v in data]

        csvout = csvutils.CsvOutput(
            self.output, self.dtypes, header=True, columns=header
        )

        csvout.write_text([','.join([str(v) for v in data])])

        # print(os.path.exists(self.output+'.yaml'))



    # =========================================================================
    # Run script
    # =========================================================================

    def main(self):

        duplication_metrics = self.extract_duplication_metrics()
        flagstat_metrics = self.extract_flagstat_metrics()
        wgs_metrics = self.extract_wgs_metrics()

        header = [
            'cell_id', 'unpaired_mapped_reads',
            'paired_mapped_reads', 'unpaired_duplicate_reads',
            'paired_duplicate_reads', 'unmapped_reads', 'percent_duplicate_reads',
            'estimated_library_size', 'total_reads', 'total_mapped_reads',
            'total_duplicate_reads', 'total_properly_paired',
            'coverage_breadth', 'coverage_depth',
        ]

        output = (self.sample_id,) + duplication_metrics + flagstat_metrics + wgs_metrics

        if self.insert_metrics:
            insert_metrics = self.extract_insert_metrics()
            output += insert_metrics
            header += ['median_insert_size',
                       'mean_insert_size',
                       'standard_deviation_insert_size']

        self.write_data(header, output)
