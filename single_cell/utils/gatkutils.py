'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import pypeliner

def generate_targets(input_bams, config, intervals, interval, **kwargs):
    # generate positions
    cmd = ['gatk', '-Xmx8G',
           '-T', 'RealignerTargetCreator',
           '-R', config['ref_genome'],
           '-o', intervals, '-L', interval,
           ]

    for _, bamfile in input_bams.items():
        cmd.extend(['-I', bamfile])

    pypeliner.commandline.execute(*cmd, **kwargs)


def gatk_realigner(inputs, config, targets, interval, tempdir, **kwargs):


    targets = os.path.abspath(targets)
    cmd = ['gatk', '-Xmx8G',
           '-T', 'IndelRealigner',
           '-R', config['ref_genome'],
           '-targetIntervals', targets,
           '--nWayOut', '_indel_realigned.bam', '-L', interval,
           '--maxReadsForRealignment','150000'
           ]

    for _, bamfile in inputs.items():
        bamfile = os.path.abspath(bamfile)
        cmd.extend(['-I', bamfile])


    cwd = os.getcwd()
    os.chdir(tempdir)

    pypeliner.commandline.execute(*cmd, **kwargs)

    os.chdir(cwd)