'''
Created on Sep 8, 2015

@author: dgrewal
'''

import os
import pandas

class MergeFiles(object):
    '''
    merges files. no overlap queries, simple concatenation
    since columns are different, select header and insert values at proper
    indices. use N/A for missing.
    '''

    def __init__(self, infile, output, typ, sep,
                 merge_type, key_cols, nan_val):
        
        self.input = infile
        self.output = output
        self.sep = sep
        self.merge_type = merge_type
        self.key_cols = key_cols
        self.type = typ
        self.nan_val = nan_val
    
    def __enter__(self):
        return self

    def __exit__(self, typ, value, traceback):
        pass

    def load(self, fname):
        '''
        load tsv file into a pandas data frame
        '''

        #return empty frame if no data
        if os.stat(fname).st_size == 0:
            return pandas.DataFrame()

        data = pandas.read_csv(fname,
                               sep=self.sep,
                               dtype={'chromosome':str, 'start': int}
                               )

        return data


    def concat_frames(self, frames):
        '''
        annotates input_df using ref_df
        '''
        out_df = pandas.concat(frames,
                               join=self.merge_type)

        return out_df


    def merge_frames(self, frames):
        '''
        annotates input_df using ref_df
        '''

        if len(frames) == 1:
            return frames[0]
        else:
            left = frames[0]
            right = frames[1]
            merged_frame = pandas.merge(left, right,
                                        how=self.merge_type,
                                        on=self.key_cols)
            for frame in frames[2:]:
                merged_frame = pandas.merge(merged_frame, frame,
                                            how=self.merge_type,
                                            on=self.key_cols)
            return merged_frame

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

        input_df.to_csv(self.output,
                        sep=self.sep,
                        index=False)


    def merge(self, frames):
        """
        merge
        """
        if self.type == 'concatenate':
            return self.concat_frames(frames)
        elif self.type == 'merge':
            return self.merge_frames(frames)


    def main(self):
        '''
        main function
        '''

        frames = []
        for fname in self.input:
            frames.append(self.load(fname))

        out_df = self.merge(frames)

        out_df = self.replace_missing_vals(out_df, nan_val=self.nan_val)

        self.write(out_df)
