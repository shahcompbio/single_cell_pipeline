'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner


def gatk_realign(inputs, outputs, targets, ref_genome, config, tempdir):


    if not os.path.exists(tempdir):
        os.makedirs(tempdir)


    new_inputs= {}    
    for key,bamfile in inputs.iteritems():
        new_bam = os.path.join(tempdir, key+'.bam')
        new_bai = os.path.join(tempdir, key+'.bam.bai')

        os.symlink(bamfile, new_bam)
        os.symlink(bamfile+'.bai', new_bai)
        new_inputs[key] = new_bam

    mapfile = os.path.join(tempdir, 'realignment_mapping.map')

    with open(mapfile, 'w') as map_outfile:
        for _,val in new_inputs.iteritems():
            val = os.path.basename(val)
            outpath = os.path.join(tempdir,val+'.realigned.bam')
            map_outfile.write(val+"\t"+outpath+"\n")

    
    cmd = [config['gatk'], '-Xmx12G',
           '-T', 'IndelRealigner',
           '-R', ref_genome,
           '-targetIntervals', targets,
           '--nWayOut', mapfile
            ]

    for _,bamfile in new_inputs.iteritems():
        cmd.extend(['-I',bamfile])
    
    pypeliner.commandline.execute(*cmd)

    with open(mapfile) as filemapping:
        for line in filemapping:
            line = line.strip().split()
            
            sample_id = line[0].strip().split('.')[0]
            bampath = line[1]
            target_path = outputs[sample_id]
            os.rename(bampath, target_path)


def realigner_target_creator(input_bams, output, ref, config):

    cmd = [config['gatk'], '-Xmx12G',
           '-T', 'RealignerTargetCreator',
           '-R', config['ref_genome'],
            '-o', output
            ]

    for _,bamfile in input_bams.iteritems():
        cmd.extend(['-I', bamfile])
    
    pypeliner.commandline.execute(*cmd)
