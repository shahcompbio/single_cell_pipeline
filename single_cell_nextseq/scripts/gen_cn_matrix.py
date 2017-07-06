'''
Created on Sep 8, 2015

@author: dgrewal
'''
import os
import argparse
import pandas as pd
import numpy as np
import warnings

class GenerateCNMatrix(object):
    '''
    merges files. no overlap queries, simple concatenation
    since columns are different, select header and insert values at proper
    indices. use N/A for missing.
    '''

    def __init__(self, args):
        self.args = args
        
        self.sep = ',' if self.args.separator == 'comma' else '\t'
    
    @staticmethod
    def replace_missing_vals(input_df, nan_val='N/A'):
        '''
        replace NaN values with nan_val
        '''
        input_df = input_df.fillna(nan_val)

        return input_df

    def write(self, input_df):
        '''
        write the dataframe to output file
        '''

        input_df.to_csv(self.args.output,
                        sep=self.sep,
                        index=False)

    def read_hmmcopy_corrected_read_file(self, sample_id):
        """
        
        """
        column_name = self.args.column_name
        data = pd.read_csv(self.args.input)
        if column_name in data.columns:
            df = data[['chr', 'start', 'end', 'width', column_name]]
        else:
            df = data[['chr', 'start', 'end', 'width']]
            
            df[column_name] = float('NaN')
            
        df = df.rename(columns = {column_name:sample_id})
        
        return df

    def read_gcbias_file(self, sample_id):
        """
        parses the gcbias data
        """
        column_name = self.args.column_name

        data = open(self.args.input).readlines()
        skiprows = [i for i,v in enumerate(data) if v[0] == '#' or v=='\n']

       #If the file is empty (only header no data) then return 0s (dummy data)
        try:
            data = pd.read_csv(self.args.input, sep='\t', skiprows=skiprows)
        except pd.io.common.EmptyDataError, e:
            warnings.warn('No data in the GCBias output')
            #If the file is empty (only header no data) then return 0s (dummy data)
            data = np.array([np.arange(100), [0]*100]).T
            data = pd.DataFrame(data, columns = ['gc', sample_id])
            return data

        data = pd.DataFrame(data[column_name])

        data['gc'] = data.index

        df = data.rename(columns={'NORMALIZED_COVERAGE':sample_id})

        df = df[['gc',sample_id]]
        return df

    def main(self):
        '''
        main function
        '''
        sample_id = self.args.sample_id
        
        if self.args.type == 'hmmcopy_corrected_reads':
            data = self.read_hmmcopy_corrected_read_file(sample_id)
        else:
            data = self.read_gcbias_file(sample_id)
        self.write(data)


def parse_args():
    '''
    specify and parse args
    '''

    parser = argparse.ArgumentParser(description='''merge tsv/csv files''')

    parser.add_argument('--input',
                        required=True,
                        help='''input files for concatenation ''')

    parser.add_argument('--sample_id',
                        required=True,
                        help='''input files for concatenation ''')

    parser.add_argument('--separator',
                            required=True,
                            default="comma",
                            choices = ("comma","tab"),
                            help='''separator type, comma for csv, tab for tsv''')

    parser.add_argument('--type',
                            required=True,
                            default="hmmcopy_corrected_reads",
                            choices = ("hmmcopy_corrected_reads","gcbias"),
                            help='''input type''')


    parser.add_argument('--column_name',
                        required=True,
                        help='''column to be used for filling the values in the matrix ''')

    parser.add_argument('--output',
                        required=True,
                        help='''path to output file''')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    ARGS = parse_args()
    m = GenerateCNMatrix(ARGS)
    m.main()


