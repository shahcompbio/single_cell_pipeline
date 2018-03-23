from __future__ import division

import argparse
import numpy as np
import pandas as pd
from copyclone.hmm import BayesianStudentsTHMM
from copyclone.mixture import BayesianMixtureOfHMMs
from statsmodels.robust.scale import stand_mad
from statsmodels.tsa.stattools import acf


def parse_args():
    """
    parses command line arguments
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('corrected_reads',
                        help='path to the gc wig file'
                        )

    parser.add_argument('reads',
                        help='path to the output csv file'
                        )

    parser.add_argument('segments',
                        help='path to the output csv file'
                        )

    parser.add_argument('metrics',
                        help='path to the output csv file'
                        )

    args = parser.parse_args()

    return args


class RunCopyClone(object):

    def __init__(self, corrected_reads, reads_out, segments, metrics):

        self.corrected_reads = corrected_reads
        self.reads_out = reads_out
        self.segments = segments
        self.metrics = metrics

        self.chromosomes = map(str, range(1, 23)) + ["X", "Y"]

    def compute_total_reads_hmmcopy(self, df):

        data = {}
        for cell, reads in df:
            total_reads_hmmcopy = sum(reads['reads'][~reads['copy'].isnull()])

            data[cell] = total_reads_hmmcopy

        return data

    def compute_mad_hmmcopy(self, df):

        data = {}
        for cell, reads in df:
            df_hmmcopy = reads[~reads['copy'].isnull()]

            if len(df_hmmcopy) > 0:
                mad_hmmcopy = stand_mad(df_hmmcopy['copy'], c=1)
            else:
                mad_hmmcopy = float('NaN')

            data[cell] = mad_hmmcopy

        return data

    def compute_mad_autosomes(self, df):

        data = {}
        for cell, reads in df:
            df_autosomes = reads[(reads['chr'] != 'X') & (reads['chr'] != 'Y')]

            if len(df_autosomes) > 0:
                mad_autosomes = stand_mad(df_autosomes['copy'], c=1)
            else:
                mad_autosomes = float('NaN')

            data[cell] = mad_autosomes
        return data

    def compute_cv_hmmcopy(self, df):

        data = {}
        for cell, reads in df:

            df_hmmcopy = reads[~reads['copy'].isnull()]

            if len(df_hmmcopy) > 0:
                df_mean = np.mean(df_hmmcopy['copy'])

                df_sd = np.std(df_hmmcopy['copy'])

                cv_hmmcopy = df_sd / df_mean

            else:
                cv_hmmcopy = float('NaN')

            data[cell] = cv_hmmcopy

        return data

    def compute_cv_neutral_state(self, df):

        data = {}
        for cell, reads in df:

            df_neutral_state = reads[
                (reads['integer_copy_number'] == 2) & (
                    ~reads['copy'].isnull())]

            if len(df_neutral_state) > 0:
                df_mean = np.mean(df_neutral_state['copy'])

                df_sd = np.std(df_neutral_state['copy'])

                cv_neutral_state = df_sd / df_mean

            else:
                cv_neutral_state = float('NaN')

            data[cell] = cv_neutral_state
        return data

    def compute_autocorrelation_hmmcopy(self, df):

        data = {}
        for cell, reads in df:

            df_hmmcopy = reads[~reads['copy'].isnull()]

            if len(df_hmmcopy) > 0:
                acf_hmmcopy = acf(df_hmmcopy['copy'], nlags=1)

                ac_hmmcopy = acf_hmmcopy[1]

            else:
                ac_hmmcopy = float('NaN')

            data[cell] = ac_hmmcopy
        return data

    def compute_mean_median_std_hmmcopy_reads_per_bin(self, df):

        mean = {}
        median = {}
        std = {}
        for cell, reads in df:

            df_hmmcopy = reads[~reads['copy'].isnull()]

            mean_hmmcopy_reads_per_bin = np.mean(df_hmmcopy['reads'])

            median_hmmcopy_reads_per_bin = np.median(df_hmmcopy['reads'])

            std_hmmcopy_reads_per_bin = np.std(df_hmmcopy['reads'])

            mean[cell] = mean_hmmcopy_reads_per_bin
            median[cell] = median_hmmcopy_reads_per_bin
            std[cell] = std_hmmcopy_reads_per_bin

        return mean, median, std

    def compute_num_empty_bins(self, df):

        data = {}
        for cell, reads in df:

            df_hmmcopy = reads[~reads['copy'].isnull()]

            empty_bins_hmmcopy = len(df_hmmcopy[df_hmmcopy['reads'] == 0])

            data[cell] = empty_bins_hmmcopy

        return data

    def compute_mad_neutral_state(self, df):
        data = {}

        for cell, reads in df:
            df_neutral_state = reads[
                (reads['integer_copy_number'] == 2) & (
                    ~reads['copy'].isnull())]

            if len(df_neutral_state) > 0:
                mad_neutral_state = stand_mad(df_neutral_state['copy'], c=1)

            else:
                mad_neutral_state = float('NaN')

            data[cell] = mad_neutral_state

        return data

    def compute_segments(self, data, state_column='integer_copy_number'):
        """
        compute segments from read data

        :param data: df with bin, copynumber
        :returns df with segments, copynumber
        """

        data = data.groupby("cell_id")

        data = {samp: data for samp, data in data}
        df_segments = pd.DataFrame()

        for sample_id in data.keys():
            if state_column in data[sample_id].columns.tolist():

                grpdata = data[sample_id].groupby(
                    'chr',
                    sort=False,
                    as_index=False)

                df_chromosomes = [
                    grpdata.get_group(chrom) for chrom in self.chromosomes]

                for df_chr in df_chromosomes:
                    index_values = df_chr.index.values.tolist()
                    change_indices = [
                        i for i in index_values[
                            1:] if df_chr[state_column][i] != df_chr[state_column][
                            i -
                            1]]

                    chrom = df_chr['chr'].unique().tolist()[0]
                    start_indices = [index_values[0]] + change_indices
                    end_indices = [
                        i - 1 for i in change_indices] + [index_values[-1]]

                    for start_index, end_index in zip(
                            start_indices, end_indices):
                        start = df_chr.ix[start_index]['start']
                        end = df_chr.ix[end_index]['end']
                        length = end - (start - 1)
                        state = df_chr.ix[
                            start_index:end_index][state_column].unique().tolist()[0]
                        median = np.nanmedian(
                            df_chr.ix[
                                start_index:end_index][state_column])

                        df_segments = df_segments.append({'chr': chrom,
                                                          'start': start,
                                                          'end': end,
                                                          'length': length,
                                                          'start_index': start_index,
                                                          'end_index': end_index,
                                                          'cell_id': sample_id,
                                                          state_column: state,
                                                          'integer_median': median}, ignore_index=True)

        return df_segments

    def fill_chromosome(self, df, column='integer_copy_number'):
        """
        apply this to single cell df where the viterbi paths have been
        mapped back to genomic coordinates (column 'copy_number'), but some
        bins are empty because they have been filtered for mappability
        :param df: dataframe with viterbipaths mapped to column
        :returns dataframe with nas filled in
        """
        df_chromosomes = [
            rows for _,
            rows in df.groupby(
                'chr',
                sort=False,
                as_index=False)]

        C = len(df_chromosomes)

        for c in range(C):
            df_chromosomes[c][column].fillna(method='ffill', inplace=True)
            df_chromosomes[c][column].fillna(method='bfill', inplace=True)

        df_concat = pd.concat(df_chromosomes, axis=0)

        return df_concat

    def initialize_naive_diploid_bayesian_hmm(self, verbose=False):
        alpha_pi = [2, 2, 50, 2, 2, 2, 2]
        alpha_A = [[1000, 2, 2, 2, 2, 2, 2],
                   [2, 1000, 2, 2, 2, 2, 2],
                   [2, 2, 10000, 2, 2, 2, 2],
                   [2, 2, 2, 1000, 2, 2, 2],
                   [2, 2, 2, 2, 1000, 2, 2],
                   [2, 2, 2, 2, 2, 1000, 2],
                   [2, 2, 2, 2, 2, 2, 1000]]
        pi = [0.05, 0.1, 0.5, 0.2, 0.05, 0.05, 0.05]
        A = [[0.994, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001],
             [0.001, 0.994, 0.001, 0.001, 0.001, 0.001, 0.001],
             [0.0001, 0.0001, 0.9994, 0.0001, 0.0001, 0.0001, 0.0001],
             [0.001, 0.001, 0.001, 0.994, 0.001, 0.001, 0.001],
             [0.001, 0.001, 0.001, 0.001, 0.994, 0.001, 0.001],
             [0.001, 0.001, 0.001, 0.001, 0.001, 0.994, 0.001],
             [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.994]]
        mu = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
        tau = [500, 25, 25, 25, 25, 25, 15]
        nu = [5, 5, 5, 5, 5, 5, 5]
        m = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
        eta = [5000, 5000, 5000, 5000, 5000, 5000, 5000]
        shape = [3, 30, 30, 30, 30, 30, 20]
        rate = [0.01, 1, 1, 1, 1, 1, 1]
        hmm_diploid = BayesianStudentsTHMM(pi, A, mu, tau, nu, alpha_pi, alpha_A, m,
                                           eta, shape, rate, name='Diploid', verbose=verbose)
        return hmm_diploid

    def initialize_naive_triploid_bayesian_hmm(self, verbose=False):
        alpha_pi = [2, 2, 2, 50, 2, 2, 2]
        alpha_A = [[1000, 2, 2, 2, 2, 2, 2],
                   [2, 1000, 2, 2, 2, 2, 2],
                   [2, 2, 1000, 2, 2, 2, 2],
                   [2, 2, 2, 10000, 2, 2, 2],
                   [2, 2, 2, 2, 1000, 2, 2],
                   [2, 2, 2, 2, 2, 1000, 2],
                   [2, 2, 2, 2, 2, 2, 1000]]
        pi = [0.05, 0.05, 0.1, 0.5, 0.2, 0.05, 0.05]
        A = [[0.994, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001],
             [0.001, 0.994, 0.001, 0.001, 0.001, 0.001, 0.001],
             [0.001, 0.001, 0.994, 0.001, 0.001, 0.001, 0.001],
             [0.0001, 0.0001, 0.0001, 0.9994, 0.0001, 0.0001, 0.0001],
             [0.001, 0.001, 0.001, 0.001, 0.994, 0.001, 0.001],
             [0.001, 0.001, 0.001, 0.001, 0.001, 0.994, 0.001],
             [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.994]]
        mu = [0.0, 1 / 3, 2 / 3, 1.0, 4 / 3, 5 / 3, 2.0]
        tau = [500, 25, 25, 25, 25, 25, 15]
        nu = [5, 5, 5, 5, 5, 5, 5]
        m = [0.0, 1 / 3, 2 / 3, 1.0, 4 / 3, 5 / 3, 2.0]
        eta = [5000, 5000, 5000, 5000, 5000, 5000, 5000]
        shape = [3, 30, 30, 30, 30, 30, 20]
        rate = [0.01, 1, 1, 1, 1, 1, 1]
        hmm_triploid = BayesianStudentsTHMM(pi, A, mu, tau, nu, alpha_pi, alpha_A, m,
                                            eta, shape, rate, name='Triploid', verbose=verbose)
        return hmm_triploid

    def initialize_naive_bayesian_mixture(self, verbose=True):
        """
        initialize a 50 50 mixture of triploid and diploid
        """
        phi = [0.5, 0.5]
        alpha_phi = [2, 2]
        hmm_diploid = self.initialize_naive_diploid_bayesian_hmm()
        hmm_triploid = self.initialize_naive_triploid_bayesian_hmm()
        hmms = [hmm_diploid, hmm_triploid]
        mixture = BayesianMixtureOfHMMs(hmms, phi, alpha_phi, verbose=verbose)
        return mixture

    def read_csv_pandas(self, infile):
        """
        read csv with cols for chr,start,end
        :params infile: csv file
        :returns pandas dataframe
        """
        data = pd.read_csv(
            infile,
            dtype={
                "chr": str,
                "start": int,
                "end": int})
        return data

    def write_csv(self, data, outfile):
        """
        write dataframe to file
        :param data: pandas dataframe
        :param outfile: file to write to.
        """
        data.to_csv(outfile, na_rep="NA", index=False)

    def get_bins_by_chromosomes(self, df):
        """
        :param pandas df with genomic coord
        :returns list of bins per chromosome, sorted
        """
        bins = []

        for chrom in self.chromosomes:
            chrom_bins = sorted(
                [v for v in df.columns.values if v[0] == chrom])
            bins.append(chrom_bins)

        return bins

    def create_dataset_copyclone(self, df, bins):
        """
        create input dataset for copyclone. must be c-contigous
        format: list of 2D matrices (one 2D matrix per chromosome)
                each matrix is cells (rows) * bins(col)
        :param df: pandas df with cells as rows and bins as cols, fill with cor_gc
        :param bins: list of bins(chrom, start, end)
        :returns list of sorted bins(chrom, start, end) per chromosome
        """
        dataset = []
        for chrom_bins in bins:
            chromdata = np.array([df[bin_name] for bin_name in chrom_bins]).T

            if not chromdata.any():
                return

            chromdata = np.ascontiguousarray(chromdata)
            dataset.append(chromdata)

        return dataset

    def parse_viterbi_paths(self, cells, viterbi_paths, bins, data):
        """
        add viterbi paths vals to the corrected data as a new col
        :param cells: list of cells, should be in same order as input dataset
        :param data: pandas df read from input corrected data,
                after adding copynumber this df gets written to file
        :param viterbi_paths: model output viterbi path matrix list(of chromosomes) of [bin x cell]
        :params bins: (chrom, start, end) in same order as copyclone dataset used for fit.
        :returns pandas dataframe with copy_number col
        """
        def get_val(data, row):
            """
            extract the copy number from viterbi path dict
            for a given row in original dataframe
            """
            if row.cell_id not in data:
                return float('nan')
            return data[row.cell_id].get(
                (row.chr, row.start, row.end), float("nan"))

        copydata = {}
        for i, chrom_bins in enumerate(bins):
            viterbi_paths_chr = viterbi_paths[i]

            assert len(cells) == len(viterbi_paths_chr)
            for cell, cell_viterbi in zip(cells, viterbi_paths_chr):

                if cell not in copydata:
                    copydata[cell] = {}

                assert len(chrom_bins) == len(cell_viterbi)
                for chrbin, vitpath in zip(chrom_bins, cell_viterbi):
                    copydata[cell][chrbin] = vitpath

        data["integer_copy_number"] = data.apply(
            lambda x: get_val(
                copydata,
                x),
            axis=1)

        return data

    def median_of_bin_residuals_from_segment_median(self, df, df_seg):
        '''
        MBRSM: This measures "dispersion", but not "integerness".
        '''
        residuals = []

        for seg_index in range(len(df_seg)):
            seg_chr = df_seg['chr'][seg_index]

            seg_start = df_seg['start'][seg_index]

            seg_end = df_seg['end'][seg_index]

            seg_median = df_seg['median'][seg_index]

            df_seg_bins = df[(df['chr'] == seg_chr) &
                             (df['start'] >= seg_start) &
                             (df['end'] <= seg_end)]

            df_seg_bins_hmmcopy = df_seg_bins[~df_seg_bins['copy'].isnull()]

            seg_residuals = [
                np.abs(
                    x -
                    seg_median) for x in df_seg_bins_hmmcopy['integer_copy_scale']]

            residuals.extend(seg_residuals)

        median_bin_residuals_median = np.median(residuals)

        return median_bin_residuals_median

    def median_of_segment_residuals_from_segment_integer(self, df_seg):
        '''
        MSRSI: This measures "integerness", but not "dispersion".
        '''
        data = {}

        for cell, segs in df_seg:

            residuals = np.abs(
                segs['integer_median'] -
                segs['integer_copy_number'])

            if len(residuals):
                median_seg_residuals_integer = np.median(residuals)
            else:
                median_seg_residuals_integer = float('nan')

            data[cell] = median_seg_residuals_integer

        return data

    def median_of_bin_residuals_from_segment_integer(self, df):
        '''
        MBRSI: This measures both "dispersion" and "integerness".
        '''
        data = {}

        for cell, reads in df:
            df_hmmcopy = reads[~reads['copy'].isnull()]

            residuals = np.abs(
                df_hmmcopy['integer_copy_scale'] -
                df_hmmcopy['integer_copy_number'])

            median_bin_residuals_integer = np.median(residuals)

            data[cell] = median_bin_residuals_integer
        return data

    def parse_labels(self, reads, segs, cells, labels, names):

        reads = reads.groupby("cell_id")

        segs = segs.groupby("cell_id")

        total_reads_hmmcopy = self.compute_total_reads_hmmcopy(reads)

        mad_hmmcopy = self.compute_mad_hmmcopy(reads)

        mad_neutral_state = self.compute_mad_neutral_state(reads)

        mad_autosomes = self.compute_mad_autosomes(reads)

        cv_hmmcopy = self.compute_cv_hmmcopy(reads)

        cv_neutral_state = self.compute_cv_neutral_state(reads)

        autocorrelation_hmmcopy = self.compute_autocorrelation_hmmcopy(reads)

        mean_hmmcopy_reads_per_bin, median_hmmcopy_reads_per_bin, std_hmmcopy_reads_per_bin = self.compute_mean_median_std_hmmcopy_reads_per_bin(
            reads)

        empty_bins_hmmcopy = self.compute_num_empty_bins(reads)

#         MBRSI_dispersion_non_integerness = self.median_of_bin_residuals_from_segment_integer(reads)

#         MBRSM_dispersion = self.median_of_bin_residuals_from_segment_median(reads, segs)
#
        MSRSI_non_integerness = self.median_of_segment_residuals_from_segment_integer(
            segs)

        metrics = pd.DataFrame()
        for i, cell in enumerate(cells):

            label = labels[i]
            name = names[i]

            cellmetrics = pd.Series({'cell_id': cell, 'label': label,
                                     'name': name,
                                     'total_reads_hmmcopy': total_reads_hmmcopy[cell],
                                     'mad_hmmcopy': mad_hmmcopy[cell],
                                     'mad_neutral_state': mad_neutral_state[cell],
                                     'mad_autosomes': mad_autosomes[cell],
                                     'cv_hmmcopy': cv_hmmcopy[cell],
                                     'cv_neutral_state': cv_neutral_state[cell],
                                     'autocorrelation_hmmcopy': autocorrelation_hmmcopy[cell],
                                     'mean_hmmcopy_reads_per_bin': mean_hmmcopy_reads_per_bin[cell],
                                     'median_hmmcopy_reads_per_bin': median_hmmcopy_reads_per_bin[cell],
                                     'std_hmmcopy_reads_per_bin': std_hmmcopy_reads_per_bin[cell],
                                     'MSRSI_non_integerness': MSRSI_non_integerness[cell],
                                     'empty_bins_hmmcopy': empty_bins_hmmcopy[cell]})

            metrics = pd.concat(
                [metrics, pd.DataFrame(cellmetrics).transpose()])

        return metrics

    def dummy_results(self, dataset):
        predicted_labels = [float('nan')] * len(dataset[0])

        predicted_names = ['NA'] * len(dataset[0])

        for chrval in dataset:
            np.place(chrval, chrval, float('nan'))

        return predicted_labels, predicted_names, dataset

    def fit(self, df):

        df_na = df.dropna(axis=1, how="any")

        if df_na.empty:
            cells = df.index
            bins = self.get_bins_by_chromosomes(df)
            dataset = self.create_dataset_copyclone(df, bins)
            predicted_labels, predicted_names, viterbi_paths = self.dummy_results(
                dataset)
        else:
            cells = df_na.index
            bins = self.get_bins_by_chromosomes(df_na)
            dataset = self.create_dataset_copyclone(df_na, bins)

            naive_mixture = self.initialize_naive_bayesian_mixture()

            naive_mixture.fit(dataset)

            predicted_labels, predicted_names, viterbi_paths = naive_mixture.predict(
                dataset)

        predicted_labels[predicted_labels == 0] = 2
        predicted_labels[predicted_labels == 1] = 3

        return cells, bins, predicted_labels, predicted_names, viterbi_paths

    def add_integer_copy_scale(self, df, cells, labels):

        df = df.groupby("cell_id")

        outdata = []
        for i, cell in enumerate(cells):
            df_cell = df.get_group(cell)
            cell_label = labels[i]

            df_cell["integer_copy_scale"] = df_cell["copy"] * cell_label

            outdata.append(df_cell)

        outdata = pd.concat(outdata)

        return outdata

    def main(self):

        data = self.read_csv_pandas(self.corrected_reads)

        # copy for running copy clone, will merge results into original later
        df = data.copy()

        df["bin"] = list(zip(df.chr, df.start, df.end))

        df = df.pivot(index='cell_id', columns='bin', values='cor_gc')

        cells, bins, predicted_labels, predicted_names, viterbi_paths = self.fit(
            df)

        data = self.parse_viterbi_paths(cells, viterbi_paths, bins, data)

        data = self.fill_chromosome(data)

        data = self.add_integer_copy_scale(data, cells, predicted_labels)

        self.write_csv(data, self.reads_out)

        segments = self.compute_segments(data)

        self.write_csv(segments, self.segments)

        metrics = self.parse_labels(
            data,
            segments,
            cells,
            predicted_labels,
            predicted_names)

        self.write_csv(metrics, self.metrics)


if __name__ == "__main__":

    args = parse_args()
    rc = RunCopyClone(
        args.corrected_reads,
        args.reads,
        args.segments,
        args.metrics)
    rc.main()
