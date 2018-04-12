'''
Created on Apr 9, 2018

@author: dgrewal
'''

import os
import shelve
import re
import fnmatch

def clean_sentinels(args):

    dirname = args["pipelinedir"]

    rundir, pattern = args["pattern"]

    rundir = os.path.join(dirname, rundir)

    if args["mode"] == "list":
        list_sentinels(rundir, pattern)
    else:
        delete_sentinels(rundir, pattern)
    
    


def list_sentinels(dirname, pattern):

    jobs_shelf = os.path.join(dirname, "jobs.shelf")
    objs_shelf = os.path.join(dirname, "objs.shelf")

    jobs = shelve.open(jobs_shelf)

    job_matches = [v for v in jobs.keys() if fnmatch.fnmatch(v, pattern)]

    objs = shelve.open(objs_shelf)

    obj_matches = [v for v in objs.keys() if fnmatch.fnmatch(v, pattern)]

    print job_matches + obj_matches
    
    
def delete_sentinels(dirname, pattern):
    
    jobs_shelf = os.path.join(dirname, "jobs.shelf")
    objs_shelf = os.path.join(dirname, "objs.shelf")

    jobs = shelve.open(jobs_shelf)

    for job in jobs.keys():
        if fnmatch.fnmatch(job, pattern):
            del jobs[job]

    jobs.close()

    objs = shelve.open(objs_shelf)

    for obj in objs.keys():
        if fnmatch.fnmatch(obj, pattern):
            del objs[obj]

    objs.close()


