'''
Created on Sep 8, 2015

@author: dgrewal
'''
import pandas as pd


class SummaryMetrics(object):
    '''
    merges files. no overlap queries, simple concatenation
    since columns are different, select header and insert values at proper
    indices. use N/A for missing.
    '''

    def __init__(self, infile, output):
        self.infile = infile
        self.output = output


    def read_csv(self, path):
        '''
        read into df
        '''
        df = pd.read_csv(path, sep=',')

        return df

    def get_mean_std(self, df, val_col, outfile, grp_col=None):
        """
        utility function to calculate and write mean and std dev of a column(val_col)
        where another col (grp_col) is equal to some value for all values
        """

        if grp_col:
            df = df.groupby(grp_col)
            for gc in df.groups.keys():
                vals = df.get_group(gc)[val_col]
                mean_ump_rds = vals.mean()
                std_ump_rds = vals.std()
    
                outfile.write(','.join([gc, str(mean_ump_rds), str(std_ump_rds)])+'\n')
        else:
            vals = df[val_col]
            mean_ump_rds = vals.mean()
            std_ump_rds = vals.std()
    
            outfile.write(','.join([str(mean_ump_rds), str(std_ump_rds)])+'\n')

    def get_cn_summary(self, df, outfile):
        '''
        '''
        #section heading and header
        outfile.write('\n#TOTAL READS MEAN AND STD DEV\n')
        outfile.write('mean,standard_deviation\n')
        self.get_mean_std(df, 'log_likelihood', outfile)

        #section heading and header
        outfile.write('\n#TOTAL READS MEAN AND STD DEV\n')
        outfile.write('mean,standard_deviation\n')
        self.get_mean_std(df, 'mad_neutral_state', outfile)

        #section heading and header
        outfile.write('\n#MSRSI_NON_INTEGERNESS MEAN AND STD DEV\n')
        outfile.write('mean,standard_deviation\n')
        self.get_mean_std(df, 'MSRSI_non_integerness', outfile)

        #section heading and header
        outfile.write('\n#MBRSM_DISPERSION MEAN AND STD DEV\n')
        outfile.write('mean,standard_deviation\n')
        self.get_mean_std(df, 'MBRSM_dispersion', outfile)

        #section heading and header
        outfile.write('\n#MBRSI_DISPERSION_NON_INTEGERNESS MEAN AND STD DEV\n')
        outfile.write('mean,standard_deviation\n')
        self.get_mean_std(df, 'MBRSI_dispersion_non_integerness', outfile)

    def get_lib_summary(self, df, outfile):
        '''
        '''

        outfile.write('#NUMBER OF LIBRARIES\n')
        outfile.write('count\n')
        outfile.write('%s\n' %str(len(df)))

        outfile.write('\n#EXPERIMENTAL_CONDITION COUNTS\n')
        outfile.write('experimental_condition,count\n')
        df_grp = df.groupby('experimental_condition')
        for ec, gp in df_grp.groups.iteritems():
            outfile.write(','.join([ec, str(len(gp))])+'\n')


        outfile.write('\n#CELL_CALL COUNTS\n')
        outfile.write('cell_call,count\n')
        df_grp = df.groupby('cell_call')
        for cc, gp in df_grp.groups.iteritems():
            outfile.write(','.join([cc, str(len(gp))])+'\n')

        

    def main(self):
        '''
        main function
        '''
        
        df = self.read_csv(self.infile)
        outfile = open(self.output, 'w')

        # library summary
        self.get_lib_summary(df, outfile)

        #copy number metric summary
        self.get_cn_summary(df, outfile)

        outfile.close()
