'''
Created on Apr 9, 2018

@author: dgrewal
'''

import fnmatch
import os

from pypeliner.sqlitedb import SqliteDb


def clean_sentinels(args):
    dirname = args["pipelinedir"]

    rundir, pattern = args["pattern"]

    rundir = os.path.join(dirname, rundir)

    if args["mode"] == "list":
        list_sentinels(rundir, pattern)
    else:
        delete_sentinels(rundir, pattern)


def list_sentinels(dirname, pattern):
    jobs_shelf = os.path.join(dirname, "jobs.db")

    jobs = SqliteDb(jobs_shelf)

    job_matches = [v for v in jobs.keys() if fnmatch.fnmatch(v, pattern)]

    jobs.close()

    matches = job_matches

    matches = '\n'.join(matches)

    print(matches)


def delete_sentinels(dirname, pattern):
    jobs_shelf = os.path.join(dirname, "jobs.db")

    jobs = SqliteDb(jobs_shelf)

    for job in jobs.keys():
        if fnmatch.fnmatch(job, pattern):
            jobs.delete(job)

    jobs.close()
