## Docker
1. You can specify the image to use by specifying `docker_image` in the context of the job. For instance: `ctx={'docker_image':'scp/single_cell_pipeline:v0.2.6'}` will use the specified image to run the task.
2. you can also specify a docker image when running a command with `pypeliner.commandline.execute` inside a job. For instance: `pypeliner.commandline.execute('bwa', 'index', 'temp/test.bam', docker_image='scp/bwa:v0.0.1')` will run the index command inside the bwa image.
3. The docker container registry server, username and password can be specified at the runtime with the `--context_config` argument. The file should be formatted as 
    ```
    docker:
      server: 'container registry server'
      username: 'insert client_id here'
      password: 'insert password here'
      mounts:
        refdata: /refdata
        datadrive: /datadrive
    ```
    The mounts must list all directories that are needed to run the pipeline.
4. The dockerfiles for the containers used by the singlecell pipeline can be found at [github](https://github.com/shahcompbio/docker_containers)

### Building containers
The pipelines can use 2 kinds of containers:
1. The containers for running the jobs submitted by pipeline. These are containers that are specified in the context section of `workflow.transform` and `workflow.commandline` calls. At the moment, these containers should have the pipeline and pypeliner installed.
2. Containers used by `pypeliner.commandline.execute` within the python function that is submitted as a job. These containers do not require any tools other than the ones needed by the job.

For instance: for the following job:

```
    workflow.commandline(
        name='extract',
        ctx={'docker_image': 'scp/test_pipeline:v0.0.1'},
        args=(
            'samtools', 'view', '-b',
            mgd.InputFile('input.bam'),
            '>',
            mgd.OutputFile('output.bam'),
        ),
    )
```
the `scp/test_pipeline:v0.0.1` container must have the pipeline code and pypeliner installed along with samtools.

On the other hand, consider a setup like the following:

```
    workflow.transform(
        name='extract',
        ctx={'docker_image': 'scp/test_pipeline:v0.0.1'},
        func=tasks.extract,
        args=(
            mgd.InputFile('input.bam'),
            mgd.OutputFile('output.bam'),
        ),
    )
```
where the python function is:
```
def extract(inputfile, outputfile):
    cmd = ['samtools', 'view', '-b',
            mgd.InputFile('input.bam'),
            '>',
            mgd.OutputFile('output.bam')]
    pypeliner.commandline.execute(*cmd, docker_image='scp/samtools:v0.0.1')
```

In this case, the `scp/test_pipeline:v0.0.1` requires pipeline code and pypeliner and `scp/samtools:v0.0.1` needs only samtools to run. 
This setup can be used to reduce the number of custom docker images. A typical pipeline with setup 2 should only require a single "base" container with the pipeline code, pypeliner and python dependencies. The third party dependencies can be run using the containers for the respective tools and can be reused across different pipelines.


### Azure Container Registry
1. Setup a docker container registry on azure. Please see [link](https://docs.microsoft.com/en-us/azure/container-registry/container-registry-get-started-portal) for details.
2. Setup an Active Directory App if you haven't already. Please see Active directory above for details.
3. Check if you can login with the active directory credentials. See [link](https://docs.microsoft.com/en-us/azure/container-registry/container-registry-authentication) for details. The `SP_APP_ID` is the CLIENT_ID and `SP_PASSWD` is the SECRET_KEY.
4. once you're logged in, you should be able to perform all operations with the docker commands.


