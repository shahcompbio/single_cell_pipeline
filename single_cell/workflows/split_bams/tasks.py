'''
Created on Nov 21, 2017

@author: dgrewal
'''
import pypeliner

def split_bam_file(bam, bai, outbam, outbai, interval):

    pypeliner.commandline.execute(
        'samtools', 'view', '-b', bam, interval,
        '>', outbam)

    pypeliner.commandline.execute(
        'samtools', 'index', outbam,
        outbai)

def split_bam_file_one_job(bam, bai, outbam, outbai, intervals):

    for interval_idx, interval in intervals.iteritems():
        output_bam = outbam(interval_idx)
        output_bai = outbai(interval_idx)

        pypeliner.commandline.execute(
            'samtools', 'view', '-b', bam, interval,
            '>', output_bam)

        pypeliner.commandline.execute(
            'samtools', 'index', output_bam,
            output_bai)