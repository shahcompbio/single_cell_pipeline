'''
Created on 2014-04-05

@author: Andrew Roth
'''
from ruffus import cmdline
from ruffus.ruffus_utility import CHECKSUM_FILE_TIMESTAMPS

def run_pipeline(args, job_manager, checksum_level=CHECKSUM_FILE_TIMESTAMPS, multithread=10, **kwargs):
    args.jobs = multithread
    
    try:
        cmdline.run(args, multithread=multithread, checksum_level=checksum_level, **kwargs)
    
    finally:
        job_manager.close()
