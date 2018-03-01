'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import pypeliner

def generate_targets(input_bams, config, intervals, interval):
    # generate positions
    cmd = ['gatk', '-Xmx8G',
           '-T', 'RealignerTargetCreator',
           '-R', config['ref_genome'],
           '-o', intervals, '-L', interval,
           ]

    for _, bamfile in input_bams.iteritems():
        cmd.extend(['-I', bamfile])

    pypeliner.commandline.execute(*cmd)


def gatk_realigner(inputs, config, targets, interval, tempdir):


    targets = os.path.abspath(targets)
    cmd = ['gatk', '-Xmx8G',
           '-T', 'IndelRealigner',
           '-R', config['ref_genome'],
           '-targetIntervals', targets,
           '--nWayOut', '_indel_realigned.bam', '-L', interval,
           '--maxReadsForRealignment','150000'
           ]

    for _, bamfile in inputs.iteritems():
        bamfile = os.path.abspath(bamfile)
        cmd.extend(['-I', bamfile])


    cwd = os.getcwd()
    os.chdir(tempdir)

    pypeliner.commandline.execute(*cmd)

    os.chdir(cwd)