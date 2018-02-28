'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import pypeliner

def generate_targets(input_bams, config, intervals, interval):
    # generate positions
    cmd = ['gatk', '-Xmx2G',
           '-T', 'RealignerTargetCreator',
           '-R', config['ref_genome'],
           '-o', intervals, '-L', interval,
           '--maxReadsForRealignment','150000'
           ]

    for _, bamfile in input_bams.iteritems():
        cmd.extend(['-I', bamfile])

    pypeliner.commandline.execute(*cmd)


def gatk_realigner(inputs, config, targets, interval, tempdir):
    cmd = ['gatk', '-Xmx2G',
           '-T', 'IndelRealigner',
           '-R', config['ref_genome'],
           '-targetIntervals', targets,
           '--nWayOut', '_indel_realigned.bam', '-L', interval,
           '--maxReadsForRealignment','150000'
           ]

    for _, bamfile in inputs.iteritems():
        cmd.extend(['-I', bamfile])

    os.chdir(tempdir)

    pypeliner.commandline.execute(*cmd)
