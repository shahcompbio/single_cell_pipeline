'''
Created on 2013-12-13

@author: Andrew Roth
'''
import shutil

def transfer_local_file(job_manager, remote_host, local_file, remote_file):
    '''
    Transfer a local file to a remote host using rsync.
    '''
    cmd = 'rsync'
    
    cmd_args = [
                '-L',
                local_file,
                '{0}:{1}'.format(remote_host, remote_file)
                ]
    
    job_manager.run_job(cmd, cmd_args)

def transfer_remote_file(job_manager, remote_host, remote_file, local_file):
    '''
    Transfer a remote file to a local host using rsync.
    '''
    tmp_file = local_file + '.tmp'
    
    cmd = 'rsync'
    
    cmd_args = [
                '-L',
                '{0}:{1}'.format(remote_host, remote_file),
                tmp_file
                ]
    
    job_manager.run_job(cmd, cmd_args)

    shutil.move(tmp_file, local_file)