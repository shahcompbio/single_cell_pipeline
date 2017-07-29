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


def merge_realignment(input_filenames, output_filename, config, input_sample_id):
    merge_filenames = []
    for (chromosome, sample_id), filename in input_filenames.iteritems():
        if input_sample_id != sample_id:
            continue
        merge_filenames.append(filename)

    merge_bams(merge_filenames, output_filename, config)


def copy_files(in_r1,out_r1):
    shutil.copy(in_r1, out_r1)


def realign(input_bams, output_bams, tempdir, config, interval):

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

    print ' '.join(cmd)
    pypeliner.commandline.execute(*cmd)

    cmd = [config['gatk'], '-Xmx12G',
           '-T', 'IndelRealigner',
           '-R', config['ref_genome'],
           '-targetIntervals', intervals,
           '--nWayOut', '_indel_realigned.bam', '-L', interval
            ]

    for _,bamfile in new_inputs.iteritems():
        cmd.extend(['-I',bamfile])
    
    os.chdir(tempdir)

    print ' '.join(cmd)
    pypeliner.commandline.execute(*cmd)

    print os.listdir(tempdir)
    for key in input_bams.keys():
        realigned_bam = os.path.join(tempdir, key+'_indel_realigned.bam')
        realigned_bai = os.path.join(tempdir, key+'_indel_realigned.bai')
        output_bam_filename = output_bams[key]
        output_bai_filename = output_bam_filename + '.bai'
        print realigned_bai, output_bai_filename
        os.rename(realigned_bam, output_bam_filename)
        os.rename(realigned_bai, output_bai_filename)



