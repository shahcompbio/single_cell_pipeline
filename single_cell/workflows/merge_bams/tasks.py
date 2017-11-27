'''
Created on Jul 24, 2017

@author: dgrewal
'''
import pypeliner

def bam_index(infile, outfile):
    pypeliner.commandline.execute(
        'samtools', 'index',
        infile,
        outfile
        )



def merge_bams(inputs, output, output_index):
    filenames = inputs.values()
    
    
    cmd = ['picard', '-Xmx2G', '-Xms2G',
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


    bam_index(output, output_index)