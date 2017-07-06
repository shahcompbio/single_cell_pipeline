'''
Created on Nov 16, 2016

@author: dgrewal
'''
from __future__ import division

import argparse
import pandas as pd
import math

def parse_args():
    #=======================================================================================================================
    # Read Command Line Input
    #=======================================================================================================================
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--corrected_reads',
                        required=True, 
                        help='''Path to HMMcopy corrected reads output .csv file.''')
    
    parser.add_argument('--segments',
                        required=True, 
                        help='''Path to HMMcopy segments output .csv file.''')
    
    parser.add_argument('--quality_metrics',
                        required=True, 
                        help='''Optional quality metrics file for the run, with 'mad_neutral_state' column.''')

    parser.add_argument('--mad_threshold', type=float, default=0,
                        help='''all cells that have low MAD won't be plotted''')

    parser.add_argument('--reads_output',
                        required=True, 
                        help='''Path to HMMcopy corrected reads output .pdf file.''')

    parser.add_argument('--segs_output',
                        required=True, 
                        help='''Path to HMMcopy segs reads output .pdf file.''')

    
    
    args = parser.parse_args()
    return args


class FilterHmmData(object):
    """
    generate the reads, bias and segment plots
    """ 
    def __init__(self, args):
        self.args = args

    def load_data_pandas(self, infile):
        """
        
        """
        data = pd.read_csv(infile,
                           sep=',')

        data = data.groupby('cell_id')
        
        return data


    def read_quality_metrics(self):
        """
        
        """
        
        df = self.load_data_pandas(self.args.quality_metrics)
        
        return df


    def read_corrected_reads(self):
        """
        
        """
        
        df = self.load_data_pandas(self.args.corrected_reads)
        
        return df
    
    def read_segments(self):
        """
        
        """
        
        df = self.load_data_pandas(self.args.segments)
        
        return df

    def get_sample_ids(self, df):
        """
        
        """
        samples = df.groups.keys()

        samples = sorted(samples)
        
        return samples

    def get_mad_score(self, sample_id, metrics):
        """
        """
        mad = metrics.get_group(sample_id)['mad_neutral_state'].iloc[0]
        return mad


    def check_mad_score(self, sample, metrics):
        """
        
        """
        mad = self.get_mad_score(sample, metrics)

        # if mad_threshold is set to nonzero.
        #zero is defaults and means mad_threshold is not set. so no filtering
        if self.args.mad_threshold:
            if math.isnan(mad):
                return False

            if mad > self.args.mad_threshold:
                return False
        return True

    def filter_write(self, df, metrics, outfile):
        """
        """
        head = False
        samples = self.get_sample_ids(df)
        
        for sample in samples:
            #If the check_mad returns false: filter it
            if not self.check_mad_score(sample, metrics):
                continue

            
            #write data
            df_samp = df.get_group(sample)
            
            if not head:
                with open(outfile, 'w') as fout:
                    df_samp.to_csv(fout, index=False)
                head = True
            else:
                with open(outfile, 'a') as fout:
                    df_samp.to_csv(fout, index=False, header=False)



    def main(self):
        """
        main
        """
        metrics = self.read_quality_metrics()

        df = self.read_corrected_reads()
        self.filter_write(df, metrics, self.args.reads_output)

        df = self.read_segments()
        self.filter_write(df, metrics, self.args.segs_output)

if __name__ == '__main__':
    args = parse_args()
    
    genhmm = FilterHmmData(args)

    genhmm.main()
