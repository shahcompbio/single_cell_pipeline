'''
Created on Jul 24, 2017

@author: dgrewal
'''
import pypeliner


def merge_bams(inputs, output):
    filenames = inputs.values()
    
    
    cmd = ['picard', '-Xmx4G',
           'MergeSamFiles',
           'OUTPUT=' + output,
           'SORT_ORDER=coordinate',
           'ASSUME_SORTED=true',
           'VALIDATION_STRINGENCY=LENIENT',
           ]
    for bamfile in filenames:
        cmd.append('I='+bamfile)
    
    pypeliner.commandline.execute(*cmd)
