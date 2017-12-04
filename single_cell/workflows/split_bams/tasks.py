'''
Created on Nov 21, 2017

@author: dgrewal
'''
import pypeliner

def split_bam_file(bam, bai, reference, outbam, outbai, interval):
    
    interval = interval.split('_')
    interval = interval[0] +':'+interval[1] +'-'+interval[2]
    
    pypeliner.commandline.execute(
        'samtools', 'view', '-b', bam, interval,
        '>', outbam)
    
    pypeliner.commandline.execute(
        'samtools', 'index', outbam,
        outbai)
