'''
Created on Jul 24, 2017

@author: dgrewal
'''
import pypeliner


def merge_bams(inputs, output):
    filenames = inputs.values()
    
    
    cmd = ['picard', '-Xmx1024m', '-Xms1024m',
           '-XX:ParallelGCThreads=1',
           'MergeSamFiles',
           'OUTPUT=' + output,
           'SORT_ORDER=coordinate',
           'ASSUME_SORTED=true',
           'VALIDATION_STRINGENCY=LENIENT',
           'MAX_RECORDS_IN_RAM=150000'
           ]
    for bamfile in filenames:
        cmd.append('I='+bamfile)
    
    pypeliner.commandline.execute(*cmd)
