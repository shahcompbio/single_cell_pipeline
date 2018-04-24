'''
Created on Apr 13, 2018

@author: dgrewal
'''
import tasks
import remixt
import pypeliner
import pypeliner.managed as mgd

default_chromosomes = [str(a) for a in xrange(1, 23)] + ['X']

def create_extract_seqdata_workflow(
     bam_filename,
     bai_filename,
     seqdata_filename,
     config,
     ref_data_dir,
):
    chromosomes = remixt.config.get_chromosomes(config, ref_data_dir)
    snp_positions_filename = remixt.config.get_filename(config, ref_data_dir, 'snp_positions')

    bam_max_fragment_length = remixt.config.get_param(config, 'bam_max_fragment_length')
    bam_max_soft_clipped = remixt.config.get_param(config, 'bam_max_soft_clipped')
    bam_check_proper_pair = remixt.config.get_param(config, 'bam_check_proper_pair')

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(obj=mgd.OutputChunks('chromosome'), value=chromosomes)

    workflow.transform(
        name='create_chromosome_seqdata',
        ctx={'mem': config["memory"]['high'], 'pool_id': config['pools']['highmem'], 'ncpus':1},
        func=tasks.create_chromosome_seqdata,
        args=(
            mgd.TempOutputFile('seqdata', 'chromosome', axes_origin=[]),
            mgd.InputFile(bam_filename),
            snp_positions_filename,
            chromosomes,
            bam_max_fragment_length,
            bam_max_soft_clipped,
            bam_check_proper_pair,
        ),
    )

    workflow.transform(
        name='merge_seqdata',
        ctx={'mem': config["memory"]['high'], 'pool_id': config['pools']['highmem'], 'ncpus':1},
        func=remixt.seqdataio.merge_seqdata,
        args=(
            mgd.OutputFile(seqdata_filename),
            mgd.TempInputFile('seqdata', 'chromosome'),
        ),
    )

    return workflow

