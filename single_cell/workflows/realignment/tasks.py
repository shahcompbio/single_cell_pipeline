'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner
import shutil


def merge_bams(inputs, output, config):

    cmd = ['picard', '-Xmx1024m', '-Xms1024m',
           '-XX:ParallelGCThreads=1',
           'MergeSamFiles',
           'OUTPUT=' + output,
           'SORT_ORDER=coordinate',
           'ASSUME_SORTED=true',
           'VALIDATION_STRINGENCY=LENIENT',
           'MAX_RECORDS_IN_RAM=1000000'
           ]
    for bamfile in inputs:
        cmd.append('I=' + bamfile)

    pypeliner.commandline.execute(*cmd)


def merge_realignment(input_filenames, output_filename,
                      config, input_sample_id):
    merge_filenames = []
    for (_, sample_id), filename in input_filenames.iteritems():
        if input_sample_id != sample_id:
            continue
        merge_filenames.append(filename)

    merge_bams(merge_filenames, output_filename, config)


def copy_files(inp, outp):
    shutil.copy(inp, outp)


def generate_targets(input_bams, config, intervals, interval):
    # generate positions
    cmd = ['gatk', '-Xmx4G',
           '-T', 'RealignerTargetCreator',
           '-R', config['ref_genome'],
           '-o', intervals, '-L', interval,
           'MAX_RECORDS_IN_RAM=1000000'
           ]

    for _, bamfile in input_bams.iteritems():
        cmd.extend(['-I', bamfile])

    pypeliner.commandline.execute(*cmd)


def gatk_realigner(inputs, config, targets, interval, tempdir):
    cmd = ['gatk', '-Xmx4G',
           '-T', 'IndelRealigner',
           '-R', config['ref_genome'],
           '-targetIntervals', targets,
           '--nWayOut', '_indel_realigned.bam', '-L', interval,
           'MAX_RECORDS_IN_RAM=1000000'
           ]

    for _, bamfile in inputs.iteritems():
        cmd.extend(['-I', bamfile])

    os.chdir(tempdir)

    pypeliner.commandline.execute(*cmd)


def realign(input_bams, input_bais, output_bams, tempdir, config, interval):

    # make the dir
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    # symlink inputs to tempdir, inputs have same filename but they should be
    # different for mapping file nwayout to work
    # realign
    new_inputs = {}
    for key, bamfile in input_bams.iteritems():
        new_bam = os.path.join(tempdir, key + '.bam')
        new_bai = os.path.join(tempdir, key + '.bam.bai')

        os.symlink(bamfile, new_bam)
        os.symlink(bamfile + '.bai', new_bai)
        new_inputs[key] = new_bam

    # save intervals file in tempdir
    targets = os.path.join(tempdir, 'realn_positions.intervals')
    generate_targets(input_bams, config, targets, interval)

    # run gatk realigner
    gatk_realigner(new_inputs, config, targets, interval, tempdir)

    # copy generated files in temp dir to the specified output paths
    for key in input_bams.keys():
        realigned_bam = os.path.join(tempdir, key + '_indel_realigned.bam')
        realigned_bai = os.path.join(tempdir, key + '_indel_realigned.bai')
        output_bam_filename = output_bams[key]
        output_bai_filename = output_bam_filename + '.bai'
        os.rename(realigned_bam, output_bam_filename)
        os.rename(realigned_bai, output_bai_filename)
