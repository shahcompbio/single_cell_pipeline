from __future__ import division
import argparse
import numpy as np
import pandas as pd
import copyclone_metrics as ccmetrics
from copyclone.hmm import BayesianStudentsTHMM
from copyclone.mixture import BayesianMixtureOfHMMs
import itertools


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

    parser.add_argument('--a',
                        default=[0.994, 0.994, 0.994, 0.994,
                                 0.994, 0.994, 0.994],
                        type=float,
                        nargs="*",
                        help='copy clone a param'
                        )

    parser.add_argument('--alpha_a',
                        default=[1000, 1000, 1000, 1000, 1000, 1000, 1000],
                        type=float,
                        nargs="*",
                        help='copy clone alpha_a param'
                        )

    parser.add_argument('--pi',
                        default=[0.05, 0.1, 0.5, 0.2, 0.05, 0.05, 0.05],
                        type=float,
                        nargs="*",
                        help='copy clone pi param'
                        )

    parser.add_argument('--alpha_pi',
                        default=[2, 2, 50, 2, 2, 2, 2],
                        type=float,
                        nargs="*",
                        help='copy clone alpha_pi param'
                        )

    parser.add_argument('--tau',
                        default=[500, 25, 25, 25, 25, 25, 15],
                        type=float,
                        nargs="*",
                        help='copy clone tau param'
                        )

    parser.add_argument('--nu',
                        default=[5, 5, 5, 5, 5, 5, 5],
                        type=float,
                        nargs="*",
                        help='copy clone nu param'
                        )

    parser.add_argument('--eta',
                        default=[5000, 5000, 5000, 5000, 5000, 5000, 5000],
                        type=float,
                        nargs="*",
                        help='copy clone params'
                        )

    parser.add_argument('--shape',
                        default=[3, 30, 30, 30, 30, 30, 20],
                        type=float,
                        nargs="*",
                        help='copy clone shape'
                        )

    parser.add_argument('--rate',
                        default=[0.01, 1, 1, 1, 1, 1, 1],
                        type=float,
                        nargs="*",
                        help='copy clone rate'
                        )

    parser.add_argument('--ploidy_states',
                        type=int,
                        nargs="*",
                        default=[2, 3, 4],
                        help='json object with all states and their corresponding labels'
                        )

    parser.add_argument('--num_states',
                        default=7,
                        type=int,
                        help='number of hmmcopy states'
                        )

    args = parser.parse_args()

    return args


class RunCopyClone(object):

    def __init__(self, corrected_reads, reads_out, segments, metrics, A=None, alpha_A=None, pi=None, alpha_pi=None,
                 tau=None, nu=None, eta=None, shape=None, rate=None, ploidy_states=None, num_states=None):

        self.corrected_reads = corrected_reads
        self.reads_out = reads_out
        self.segments = segments
        self.metrics = metrics

        self.A = A
        self.alpha_A = alpha_A
        self.pi = pi
        self.alpha_pi = alpha_pi
        self.tau = tau
        self.nu = nu
        self.eta = eta
        self.shape = shape
        self.rate = rate
        self.ploidy_states = ploidy_states
        self.num_states = num_states

        self.chromosomes = map(str, range(1, 23)) + ["X", "Y"]

    def compute_segments(self, data, state_column='state'):
        """
        compute segments from read data

        :param data: df with bin, copynumber
        :returns df with segments, copynumber
        """
        df_segments = pd.DataFrame()

        assert state_column in data

        data = data.groupby(["cell_id", "chr"])

        for (cellid, chrom), chromdata in data:

            chromdata = chromdata.reset_index()

            i = 0
            for _, label in itertools.groupby(chromdata[state_column]):
                label = list(label)
                start_index = i
                end_index = i + len(label) - 1
                i += len(label)

                start = chromdata.at[start_index, 'start']
                end = chromdata.at[end_index, 'end']
                length = end - (start - 1)
                state = label[0]

                # assuming copy is already scaled by ploidy
                median = np.nanmedian(chromdata[start_index:end_index]['copy'])

                df_segments = df_segments.append({'chr': chrom,
                                                  'start': start,
                                                  'end': end,
                                                  'length': length,
                                                  'start_index': start_index,
                                                  'end_index': end_index,
                                                  'cell_id': cellid,
                                                  state_column: state,
                                                  'median': median}, ignore_index=True)

        return df_segments

    def fill_chromosome(self, df, column='state'):
        """
        apply this to single cell df where the viterbi paths have been
        mapped back to genomic coordinates (column 'copy_number'), but some
        bins are empty because they have been filtered for mappability
        :param df: dataframe with viterbipaths mapped to column
        :returns dataframe with nas filled in
        """
        df_chromosomes = [
            rows for _, rows in df.groupby('chr', sort=False, as_index=False)
        ]

        for c in range(len(df_chromosomes)):
            df_chromosomes[c][column].fillna(method='ffill', inplace=True)
            df_chromosomes[c][column].fillna(method='bfill', inplace=True)

        df_concat = pd.concat(df_chromosomes, axis=0)

        return df_concat

    def get_param_matrix(self, vals, num_states, fill_val=None):

        if not fill_val:
            param = [
                [(1 - val) / (num_states - 1)] * num_states for val in vals]
            param = np.matrix(param)
        else:
            param = np.full([num_states, num_states], fill_val)

        np.fill_diagonal(param, vals)

        return param

    def initialize_naive_bayesian_hmm(self, verbose=False):

        A = self.get_param_matrix(self.A, self.num_states)
        alpha_A = self.get_param_matrix(
            self.alpha_A,
            self.num_states,
            fill_val=2)

        hmms = []

        for val in self.ploidy_states:

            mu = [0.0, 1 / val, 2 / val, 3 / val, 4 / val, 5 / val, 6 / val]
            m = [0.0, 1 / val, 2 / val, 3 / val, 4 / val, 5 / val, 6 / val]

            print (self.pi, A, mu, self.tau, self.nu, self.alpha_pi, alpha_A, m,
                   self.eta, self.shape, self.rate)

            hmm = BayesianStudentsTHMM(
                self.pi, A, mu, self.tau, self.nu, self.alpha_pi, alpha_A, m,
                self.eta, self.shape, self.rate, name=str(val), verbose=verbose
            )

            hmms.append(hmm)

        return hmms

    def initialize_naive_bayesian_mixture(self, verbose=True):
        """
        initialize a 50 50 mixture of triploid and diploid
        """

        num_states = len(self.ploidy_states)

        phi = [1 / num_states] * num_states
        alpha_phi = [2] * num_states

        hmms = self.initialize_naive_bayesian_hmm(verbose=verbose)
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

    def parse_viterbi_paths(
            self, cells, viterbi_paths, bins, data, state_col="state"):
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

        data[state_col] = data.apply(lambda x: get_val(copydata, x), axis=1)

        return data

    def parse_labels(self, reads, segs, cells, names):

        reads = reads.groupby("cell_id")

        segs = segs.groupby("cell_id")

        total_reads_hmmcopy = ccmetrics.compute_total_reads_hmmcopy(reads)

        mad_hmmcopy = ccmetrics.compute_mad_hmmcopy(reads)

        mad_neutral_state = ccmetrics.compute_mad_neutral_state(reads)

        mad_autosomes = ccmetrics.compute_mad_autosomes(reads)

        cv_hmmcopy = ccmetrics.compute_cv_hmmcopy(reads)

        cv_neutral_state = ccmetrics.compute_cv_neutral_state(reads)

        autocorrelation_hmmcopy = ccmetrics.compute_autocorrelation_hmmcopy(
            reads)

        mean_hmmcopy_reads_per_bin, median_hmmcopy_reads_per_bin, std_hmmcopy_reads_per_bin = ccmetrics.compute_mean_median_std_hmmcopy_reads_per_bin(
            reads)

        empty_bins_hmmcopy = ccmetrics.compute_num_empty_bins(reads)

        MBRSI_dispersion_non_integerness = ccmetrics.median_of_bin_residuals_from_segment_integer(
            reads)

        MBRSM_dispersion = ccmetrics.median_of_bin_residuals_from_segment_median(
            reads,
            segs)

        MSRSI_non_integerness = ccmetrics.median_of_segment_residuals_from_segment_integer(
            segs)

        halfiness = ccmetrics.compute_halfiness(reads, segs)

        metrics = pd.DataFrame()
        for i, cell in enumerate(cells):

            name = names[i]

            cellmetrics = pd.Series({'cell_id': cell,
                                     'ploidy': name,
                                     'total_mapped_reads': total_reads_hmmcopy[cell],
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
                                     'MBRSI_dispersion_non_integerness': MBRSI_dispersion_non_integerness[cell],
                                     'MBRSM_dispersion': MBRSM_dispersion[cell],
                                     'empty_bins_hmmcopy': empty_bins_hmmcopy[cell],
                                     'total_halfiness': halfiness[cell][0],
                                     'scaled_halfiness': halfiness[cell][1],
                                     })

            metrics = pd.concat(
                [metrics, pd.DataFrame(cellmetrics).transpose()])

        metrics["log_likelihood"] = float("nan")

        return metrics

    def dummy_results(self, dataset):
        predicted_labels = [float('nan')] * len(dataset[0])

        predicted_names = ['nan'] * len(dataset[0])

        for chrval in dataset:
            np.place(chrval, chrval, float('nan'))

        return predicted_labels, predicted_names, dataset

    def scale_copy_by_ploidy(self, df, cells, labels):

        df = df.groupby("cell_id")

        outdata = []
        for i, cell in enumerate(cells):
            df_cell = df.get_group(cell)
            cell_label = labels[i]

            df_cell["copy"] = df_cell["copy"] * cell_label

            outdata.append(df_cell)

        outdata = pd.concat(outdata)

        return outdata

    def fit(self, df):

        df_na = df.dropna(axis=1, how="any")

        if df_na.empty:
            cells = df.index
            bins = self.get_bins_by_chromosomes(df)
            dataset = self.create_dataset_copyclone(df, bins)
            _, predicted_names, viterbi_paths = self.dummy_results(
                dataset)
        else:
            cells = df_na.index
            bins = self.get_bins_by_chromosomes(df_na)

            dataset = self.create_dataset_copyclone(df_na, bins)

            naive_mixture = self.initialize_naive_bayesian_mixture()

            naive_mixture.fit(dataset)

            _, predicted_names, viterbi_paths = naive_mixture.predict(
                dataset)

        predicted_names = map(float, predicted_names)

        return cells, bins, predicted_names, viterbi_paths

    def main(self):

        data = self.read_csv_pandas(self.corrected_reads)

        # copy for running copy clone, will merge results into original later
        df = data.copy()

        df["bin"] = list(zip(df.chr, df.start, df.end))

        df = df.pivot(index='cell_id', columns='bin', values='cor_gc')

        cells, bins, predicted_names, viterbi_paths = self.fit(
            df)

        data = self.parse_viterbi_paths(cells, viterbi_paths, bins, data)

        data = self.fill_chromosome(data)

        data = self.scale_copy_by_ploidy(data, cells, predicted_names)

        self.write_csv(data, self.reads_out)

        segments = self.compute_segments(data)

        self.write_csv(segments, self.segments)

        metrics = self.parse_labels(
            data,
            segments,
            cells,
            predicted_names)

        self.write_csv(metrics, self.metrics)


if __name__ == "__main__":

    args = parse_args()

    rc = RunCopyClone(args.corrected_reads, args.reads, args.segments, args.metrics,
                      A=args.a, alpha_A=args.alpha_a, pi=args.pi, alpha_pi=args.alpha_pi,
                      tau=args.tau, nu=args.nu, eta=args.eta, shape=args.shape, rate=args.rate,
                      ploidy_states=args.ploidy_states, num_states=args.num_states)
    rc.main()
