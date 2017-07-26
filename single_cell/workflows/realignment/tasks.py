'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner
import shutil

def merge_bams(inputs, output, config):
    
    cmd = [config['picard'], '-Xmx12G',
           'MergeSamFiles',
           'OUTPUT=' + output,
           'SORT_ORDER=coordinate',
           'ASSUME_SORTED=true',
           'VALIDATION_STRINGENCY=LENIENT',
           ]
    for bamfile in inputs:
        cmd.append('I='+bamfile)
    
    pypeliner.commandline.execute(*cmd)



def merge_realignment(indir, output, config, sample_id):

    bams = [os.path.join(v, sample_id+'.bam') for v in indir.values()]

    merge_bams(bams, output, config)



def copy_files(in_r1,out_r1):
    shutil.copy(in_r1, out_r1)

def realign(input_bams, tempdir, config, interval):

    #make the dir
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    #symlink inputs to tempdir, inputs have same filename but they should be
    #different for mapping file nwayout to work
    #realign
    new_inputs= {}    
    for key,bamfile in input_bams.iteritems():
        new_bam = os.path.join(tempdir, key+'.bam')
        new_bai = os.path.join(tempdir, key+'.bam.bai')

        os.symlink(bamfile, new_bam)
        os.symlink(bamfile+'.bai', new_bai)
        new_inputs[key] = new_bam

    #create mapping file for realign output
    mapfile = os.path.join(tempdir, 'realignment_mapping.map')

    with open(mapfile, 'w') as map_outfile:
        for _,val in new_inputs.iteritems():
            val = os.path.basename(val)
            outpath = os.path.join(tempdir,val+'.realigned.bam')
            map_outfile.write(val+"\t"+outpath+"\n")

    #save intervals file in tempdir
    intervals = os.path.join(tempdir, 'realn_positions.intervals')

    #generate positions    
    cmd = [config['gatk'], '-Xmx12G',
           '-T', 'RealignerTargetCreator',
           '-R', config['ref_genome'],
            '-o', intervals, '-L', interval
            ]

    for _,bamfile in input_bams.iteritems():
        cmd.extend(['-I', bamfile])
    
    pypeliner.commandline.execute(*cmd)

    cmd = [config['gatk'], '-Xmx12G',
           '-T', 'IndelRealigner',
           '-R', config['ref_genome'],
           '-targetIntervals', intervals,
           '--nWayOut', mapfile, '-L', interval
            ]

    for _,bamfile in new_inputs.iteritems():
        cmd.extend(['-I',bamfile])
    
    pypeliner.commandline.execute(*cmd)

