from pipelines.io import make_directory

import os
import subprocess
import time

class JobManager(object):
    def run_job(self, cmd, cmd_args, **kwargs):
        '''
        Run a job and wait for it to complete.
        '''
        pass    
 
    def close(self):
        '''
        Close the job manager.
        '''
        pass
    
class LocalJobManager(JobManager):
    def __init__(self):
        self.cmd_strs = []
 
    def run_job(self, cmd, cmd_args, **kwargs):
        cmd_str = ' '.join([cmd, ] + [str(x) for x in cmd_args])
        
        self._run_cmd(cmd_str)
            
    def _run_cmd(self, cmd_str):
        '''
        Throw exception if run command fails
        '''
        process = subprocess.Popen(cmd_str, stdin=subprocess.PIPE, shell=True)
       
        process.stdin.close()
        
        sts = os.waitpid(process.pid, 0)    
        
        if sts[1] != 0:
            raise Exception('Failed to run {0}\n'.format(cmd_str))

class ClusterJobManager(JobManager):
    def __init__(self, log_dir=None, max_tries=1):
        import drmaa
        
        self._drmaa = drmaa
        
        self.log_dir = log_dir
        
        self.max_tries = max_tries
        
        if self.log_dir is not None:
            make_directory(self.log_dir)
        
        self._init_session()

    def _init_session(self):
        self.session = self._drmaa.Session()
        
        session_started = False
        
        max_tries = 10
        
        num_tries = 0
        
        while not session_started:
            num_tries += 1
            
            try:
                self.session.initialize()
                
                session_started = True
            
            except self._drmaa.errors.DrmCommunicationException:
                if num_tries >= max_tries:
                    raise
                
                else:
                    time.sleep(60)
   
    def run_job(self, cmd, cmd_args, cleanup_log_files=False, job_name=None, mem=1, max_mem=1, num_cpus=1, queue='all.q', hosts=None):
        job = ClusterJob(cmd, cmd_args, mem, max_mem, num_cpus, queue, job_name=job_name, hosts=hosts)
        
        num_tries = 0
        
        while True:
            num_tries += 1
            
            job_id = self._run_job(job)
            
            err_file, log_file = self._get_output_files(job, job_id)
            
            try:
                self._check_exit_status(job_id)
                
                if cleanup_log_files:
                    os.unlink(err_file)
                    
                    os.unlink(log_file)
                
                break
            
            except DrmaaCommunicationException:
                self.session.control(job_id, self._drmaa.JobControlAction.TERMINATE)
            
            except Exception as e:
                if num_tries >= self.max_tries:
                    cmd_str = ' '.join([cmd, ] + [str(x) for x in cmd_args])
                    
                    exception_str = (
                                     e.message,
                                     '#' * 120,
                                     'cmd: {0}'.format(cmd_str),
                                     'stderr_file: {0}'.format(err_file),
                                     'stdout_file: {0}'.format(log_file),
                                     '#' * 120
                                     )
                    
                    exception_str = '\n'.join(exception_str)
                    
                    raise Exception(exception_str)
                
                elif cleanup_log_files:
                    os.unlink(err_file)
                    
                    os.unlink(log_file)
                
    
    def close(self):
        self.session.control(self._drmaa.Session.JOB_IDS_SESSION_ALL, self._drmaa.JobControlAction.TERMINATE)
        
        self.session.exit()

    def _check_exit_status(self, job_id):        
        job_info = self.session.wait(job_id, self._drmaa.Session.TIMEOUT_WAIT_FOREVER)
        
        if job_info.wasAborted:
            raise ClusterException('Job {0} was aborted.'.format(job_info.jobId))
        
        elif job_info.hasExited:
            job_str = 'Job {0} finished regularly with exit status {1}.'.format(job_info.jobId, job_info.exitStatus)
            
            if job_info.exitStatus == 0:
                print 'Job {0} finished regularly with exit status {1}.'.format(job_info.jobId, job_info.exitStatus)
            
            else:
                raise JobException(job_str)
        
        elif job_info.hasSignal:
            exception_str = 'Job {0} finished due to signal {1}.'.format(job_info.jobId, job_info.terminatedSignal)
        
            if job_info.hasCoreDump:
                exception_str += '\n'
                
                exception_str += 'A core dump is available for job {0}.'.format(job_info.jobId)
            
            raise ClusterException(exception_str)
         
        else:
            raise ClusterException('Job {0} finished with unclear conditions.'.format(job_info.jobId))
        
        self._print_resource_usage(job_info)
         
        return job_info.exitStatus
        
    def _run_job(self, job):
        job_submitted = False
        
        while not job_submitted:
            job_template = job.get_job_template(self.log_dir, self.session)
            try:
                job_id = self.session.runJob(job_template)
                
                job_submitted = True
                
            except self._drmaa.errors.DrmCommunicationException:
                time.sleep(60)
                
            finally:
                self.session.deleteJobTemplate(job_template)
        
        return job_id
    
    def _print_resource_usage(self, job_info):
        print '-' * 120
        print '# Job {0} resource usage.'.format(job_info.jobId)
        print '-' * 120
        
        for key, value in sorted(job_info.resourceUsage.items()):
            print key, value
    
    def _get_output_files(self, job, job_id):
        if self.log_dir is None:
            log_dir = os.getcwd()
            
        else:
            log_dir = self.log_dir
        
        if job.job_name is not None:
            prefix = job.job_name
        
        else:
            prefix = job.cmd
        
        err_file = os.path.join(log_dir, '{0}.e{1}'.format(prefix, job_id))
        
        log_file = os.path.join(log_dir, '{0}.o{1}'.format(prefix, job_id))
        
        return err_file, log_file
        
class ClusterJob(object):
    def __init__(self, cmd, cmd_args, mem, max_mem, num_cpus, queue, job_name=None, hosts=None):
        self.cmd = cmd
        
        self.cmd_args = cmd_args
        
        self.mem = mem
        
        self.max_mem = max_mem
        
        self.num_cpus = num_cpus
        
        self.job_name = job_name

        self.queue = queue 

        self.hosts = hosts

    def get_job_template(self, log_dir, session):
        job_template = session.createJobTemplate()
        
        job_template.remoteCommand = self.cmd
        
        job_template.args = [str(x) for x in self.cmd_args]
        
        job_template.workingDirectory = os.getcwd()
        
        if self.job_name is not None:
            job_template.jobName = self.job_name
        
        if log_dir is not None:
            job_template.errorPath = ':' + log_dir
            
            job_template.outputPath = ':' + log_dir
        
        if self.num_cpus > 1:
            if self.hosts:
                job_template.nativeSpecification = ' -hard -q {queue} -l mem_free={mem}G,mem_token={mem}G,h_vmem={max_mem}G -V -w n -pe ncpus {num_cpus} -l {hosts}'.format(**self.__dict__)
            else:
                job_template.nativeSpecification = ' -hard -q {queue} -l mem_free={mem}G,mem_token={mem}G,h_vmem={max_mem}G -V -w n -pe ncpus {num_cpus}'.format(**self.__dict__)
        
        else:
            if self.hosts:
                job_template.nativeSpecification = ' -hard -q {queue} -l mem_free={mem}G,mem_token={mem}G,h_vmem={max_mem}G -V -w n -l {hosts}'.format(**self.__dict__)
            else:
                job_template.nativeSpecification = ' -hard -q {queue} -l mem_free={mem}G,mem_token={mem}G,h_vmem={max_mem}G -V -w n'.format(**self.__dict__)

        return job_template

class RemoteJobManager(object):
    def __init__(self, host_name, ssh_private_key, ssh_user_name):
        import spur
        
        self._spur = spur
        
        self.host_name = host_name
        
        self.ssh_private_key = ssh_private_key
        
        self.ssh_user_name = ssh_user_name
        
        self.cmds = []
        
    def add_job(self, cmd, cmd_args, **kwargs):
        self.cmds.append((cmd, cmd_args))
    
    def run_job(self, cmd, cmd_args, **kwargs):
        self._run_cmd(cmd, cmd_args)
    
    def wait(self):
        for cmd, cmd_args in self.cmds:
            self._run_cmd(cmd, cmd_args)        
        
    def _run_cmd(self, cmd, cmd_args):
        shell = self._get_shell()
        
        results = shell.run([cmd, ] + cmd_args)
        
        if results.return_code != 0:
            raise results.to_error()

    def _get_shell(self):
        shell = self._spur.SshShell(hostname=self.host_name,
                                    private_key_file=self.ssh_private_key,
                                    username=self.ssh_user_name)
        
        return shell
    
class ClusterException(Exception):
    pass

class JobException(Exception):
    pass

class DrmaaCommunicationException(Exception):
    pass

#=======================================================================================================================
# Helper functions
#=======================================================================================================================
def get_job_managers(run_local, log_dir=None, max_tries=3):
    '''
    Initialise a session with a local job manager and a cluster job manager.
    '''
    local_job_manager = LocalJobManager()
    
    run_local_cmd = local_job_manager.run_job
    
    if run_local:
        job_manager = local_job_manager
    
    else:
        job_manager = ClusterJobManager(log_dir=log_dir, max_tries=max_tries)
    
    run_cmd = job_manager.run_job
    
    return run_local_cmd, run_cmd, local_job_manager, job_manager
