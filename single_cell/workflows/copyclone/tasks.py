'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner
from scripts import ReadCounter
from scripts import CorrectReadCount
from scripts import RunCopyClone


from single_cell.utils import csvutils
from single_cell.utils import pdfutils
from single_cell.utils import helpers

scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
run_hmmcopy_rscript = os.path.join(scripts_directory, 'hmmcopy.R')

def merge_reads(reads, merged_reads):
    csvutils.concatenate_csv_lowmem(reads, merged_reads)

def correct_reads(
    bam_file,
    bai_file,
    reads_filename,
    config,
    tempdir,
    sample_id):


    ccparams = config["copyclone"]
    run_readcount_rscript = os.path.join(scripts_directory, 'correct_read_count.R')

    #generate wig file for hmmcopy
    os.makedirs(tempdir)
    readcount_wig = os.path.join(tempdir, 'readcounter.wig')

    rc = ReadCounter(bam_file, readcount_wig, ccparams['bin_size'], config['chromosomes'],
                     ccparams['min_mqual'], excluded=ccparams['exclude_list'])
    rc.main()

    if ccparams["smoothing_function"] == 'loess':
        cmd=['Rscript', run_readcount_rscript,
             readcount_wig,
             ccparams['gc_wig_file'],
             ccparams['map_wig_file'],
             reads_filename,
             sample_id
             ]
        pypeliner.commandline.execute(*cmd)
    elif ccparams["smoothing_function"] == 'modal':
        CorrectReadCount(ccparams["gc_wig_file"],
                         ccparams['map_wig_file'],
                         readcount_wig,
                         reads_filename,
                         mappability=ccparams['map_cutoff'],
                         sample_id = sample_id).main()
    else:
        raise Exception("smoothing function %s not supported. pipeline supports loess and modal" %ccparams["smoothing_function"])


def run_copyclone(corrected_data, reads, segments, metrics):
    rc = RunCopyClone(corrected_data, reads, segments, metrics)
    rc.main()

