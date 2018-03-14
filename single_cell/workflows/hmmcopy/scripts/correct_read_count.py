'''
Created on Feb 21, 2018

@author: dgrewal
'''
from __future__ import division

import argparse
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats.mstats import mquantiles


class CorrectReadCount(object):
    """
    fit lowess/polynomial curve smoothing to reads-gc
    use fitted model to predict corrected gc and mappability
    values
    """

    def __init__(self, gc, mapp, wig, output, mappability=0.9,
                 smoothing_function='lowess',
                 polynomial_degree=2):
        self.mappability = mappability

        self.gc = gc
        self.mapp = mapp
        self.wig = wig
        self.output = output

    def read_wig(self, infile, counts=False):
        """read wiggle files

        :param infile: input wiggle file
        :param counts: set to true if infile wiggle has integer values
        """

        data = []

        with open(infile) as wig:
            for line in wig:
                line = line.strip()

                if line.startswith('fixedStep'):
                    line = line.strip().split()

                    chrom = line[1].split('=')[1]
                    winsize = int(line[3].split('=')[1])
                    start = int(line[2].split('=')[1])

                    bin_start = 0 if start < winsize else start / winsize
                else:
                    value = int(line) if counts else float(line)
                    data.append((chrom, (bin_start * winsize) + 1,
                                 (bin_start + 1) * winsize, winsize, value))
                    bin_start += 1

        return data

    def valid(self, df):
        """adds valid column (calls with atleast one reads and non negative gc)

        :params df: pandas dataframe
        """

        df.loc[:, "valid"] = True

        df.loc[(df["reads"] <= 0) | (df['gc'] < 0), "valid"] = False

        return df

    def ideal(self, df):
        """adds ideal column

        :params df: pandas dataframe
        """
        df.loc[:, "ideal"] = True

        valid_reads = df[df["valid"]]["reads"]
        valid_gc = df[df["valid"]]["gc"]

        routlier = 0.01
        doutlier = 0.001

        range_l, range_h = mquantiles(valid_reads, prob=[0, 1 - routlier],
                                      alphap=1, betap=1)
        domain_l, domain_h = mquantiles(valid_gc, prob=[doutlier, 1 - doutlier],
                                        alphap=1, betap=1)

        df.loc[(df["valid"] == False) |
               (df["map"] < self.mappability) |
               (df["reads"] <= range_l) |
               (df["reads"] > range_h) |
               (df["gc"] < domain_l) |
               (df["gc"] > domain_h),
               "ideal"] = False

        return df


    def create_dataframe(self, reads, mapp, gc):
        """merge data from reads, mappability and gc wig files
        into pandas dataframe

        :param reads: list of tuples, formatted as [(chromosome,
                    start, end, count), ]
        :param mapp: list of tuples, formatted as [(chromosome,
                    start, end, mappability value), ]
        :param reads: list of tuples, formatted as [(chromosome,
                    start, end, gc content), ]
        """
        err_str = 'please ensure that reads, mappability and '\
            'gc wig files have the same sort order'

        data = []
        for read_v, mapp_v, gc_v in zip(reads, mapp, gc):
            assert read_v[0] == mapp_v[0] == gc_v[0], err_str
            assert read_v[1] == mapp_v[1] == gc_v[1], err_str
            assert read_v[2] == mapp_v[2] == gc_v[2], err_str
            assert read_v[3] == mapp_v[3] == gc_v[3], err_str

            data.append((read_v[0], read_v[1], read_v[2], read_v[3], gc_v[4],
                         mapp_v[4], read_v[4],))

        labels = ['chr', 'start', 'end', 'width', 'gc', 'map', 'reads']
        data = pd.DataFrame(data, columns=labels)

        return data


    def modal_quantile_regression(self, df_regression, lowess_frac=0.2):
        '''
        Compute quantile regression curves and select the modal quantile.
        '''
        # 2nd order polynomial quantile regression 
        
        q_range = range(10, 91, 1)
        quantiles = np.array(q_range) / 100
        quantile_names = [str(x) for x in q_range]

        #need at least 3 values to compute the quantiles
        if len(df_regression) < 10:
            return df_regression

        poly2_quantile_model = smf.quantreg('reads ~ gc + I(gc ** 2.0)', data=df_regression)
        poly2_quantile_fit = [poly2_quantile_model.fit(q=q) for q in quantiles]
        poly2_quantile_predict = [poly2_quantile_fit[i].predict(df_regression) for i in range(len(quantiles))]
        
        poly2_quantile_params = pd.DataFrame()
        
        for i in range(len(quantiles)):
            df_regression[quantile_names[i]] = poly2_quantile_predict[i]
            poly2_quantile_params[quantile_names[i]] = poly2_quantile_fit[i].params
        
        # integration and mode selection
        
        gc_min = df_regression['gc'].quantile(q=0.10)
        gc_max = df_regression['gc'].quantile(q=0.90)
        
        poly2_quantile_integration = np.zeros(len(quantiles)+1)
        
        for i in range(len(quantiles)):
            params = poly2_quantile_params[quantile_names[i]].tolist()
            params.reverse()
            poly2 = np.poly1d(params)
            integ = poly2.integ()
            integrand = integ(gc_max) - integ(gc_min)
            poly2_quantile_integration[i+1] = integrand
        
        # find the modal quantile
        
        distances = poly2_quantile_integration[1:] - poly2_quantile_integration[:-1]
        
        df_dist = pd.DataFrame({'quantiles': quantiles, 'quantile_names': quantile_names, 'distances': distances})
        dist_max = df_dist['distances'].quantile(q=0.95)
        df_dist_filter = df_dist[df_dist['distances']<dist_max]
        df_dist_filter['lowess'] = lowess(df_dist_filter['distances'], df_dist_filter['quantiles'], frac=lowess_frac, return_sorted=False)
        
        modal_quantile = quantile_names[np.argmin(df_dist_filter['lowess'])]
        
        # add values to table
        
        df_regression['modal_quantile'] = modal_quantile
        df_regression['modal_curve'] = df_regression[modal_quantile]
        df_regression['modal_corrected'] = df_regression['reads'] / df_regression[modal_quantile]
        
        return df_regression

    def write(self, df):
        """write results to the output file

        :param df: pandas dataframe
        """

        df.to_csv(self.output, index=False, sep=',', na_rep="NA")

    def main(self):
        gc = self.read_wig(self.gc)
        mapp = self.read_wig(self.mapp)
        reads = self.read_wig(self.wig, counts=True)

        df = self.create_dataframe(reads, mapp, gc)

        df = self.valid(df)
        df = self.ideal(df)

        df['modal_quantile'] = 'NaN'
        df['modal_curve'] = 'NaN'
        df['modal_corrected'] = 'NaN'

        # filtering and sorting
        df_valid_gc = df[df['gc']>0]

        df_non_zero = df_valid_gc[df_valid_gc['reads']>0]

        df_regression = pd.DataFrame.copy(df_non_zero)

        df_regression.sort_values(by='gc', inplace=True)

        # modal quantile regression
        df_regression = self.modal_quantile_regression(df_regression, lowess_frac=0.2)

        # map results back to full data frame
        df.ix[df_regression.index, 'modal_quantile'] = df_regression['modal_quantile']
        df.ix[df_regression.index, 'modal_curve'] = df_regression['modal_curve']
        df.ix[df_regression.index, 'modal_corrected'] = df_regression['modal_corrected']
        
        # filter by mappability
        df['copy'] = df['modal_corrected']
        df['copy'][df['map'] < self.mappability] = float('NaN')

        df = df.rename(columns=({ "modal_corrected" : "cor_gc"}))

        df["cor_map"] = float("NaN")

        # save
        self.write(df)




def parse_args():
    """
    parses command line arguments
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('gc',
                        help='path to the gc wig file'
                        )

    parser.add_argument('map',
                        help='path to the mappability wig file'
                        )

    parser.add_argument('reads',
                        help='path to the read-counts wig file'
                        )

    parser.add_argument('output',
                        help='path to the output csv file'
                        )

    parser.add_argument('--mappability',
                        default=0.9,
                        type=float,
                        help='specify mappability threshold')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = parse_args()

    corr = CorrectReadCount(args.gc, args.map, args.reads, args.output,
                            mappability=args.mappability,
                            )

    corr.main()
