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

def split_bam_file_one_job(bam, bai, reference, outbam, outbai, intervals):

    for interval in intervals:
        output_bam = outbam(interval)
        output_bai = outbai(interval)

        interval = interval.split('_')
        interval = interval[0] +':'+interval[1] +'-'+interval[2]

        pypeliner.commandline.execute(
            'samtools', 'view', '-b', bam, interval,
            '>', output_bam)

        pypeliner.commandline.execute(
            'samtools', 'index', output_bam,
            output_bai)