
## Azure Batch
This guide assumes that the input files have been uploaded to azure blob storage and assumes that the jobs will be run with docker. Please see the respective guides to set them up.
1. Create a new batch account in azure. Please follow `create-a-batch-account` section at [link](https://docs.microsoft.com/en-us/azure/batch/quick-create-portal#create-a-batch-account) for details.
2. Build a Managed Image with the required dependencies. see below
3. Start a small ubuntu VM from the captured image in step 2, this machine will serve as our master node. Please see [link](https://docs.microsoft.com/en-us/azure/virtual-machines/windows/create-vm-generalized-managed) for details.
4. log into the azure container registry. see Docker section for details
5. run the single cell pipeline with docker. See below


### Launch the pipeline with docker
##### run the pipeline inside a docker container:
```
docker run shahlab.azurecr.io/scp/single_cell_pipeline:v0.2.6 single_cell -h
```
Dont forget to mount the directories and set the working directory!
##### Launch the pipeline with `--run_with_docker` flag.
There is another way to launch the pipeline with docker that requires some extra steps. 
1. Install miniconda.
2. install single_cell_pipeline and pypeliner.
3. run the single_cell_pipeline as you normally would along with the `--run_with_docker` flag.

### Build Azure Managed Image

1. start a ubuntu 18.04 LTS. 
2. Install docker on the machine. Please see [link](https://docs.docker.com/install/linux/docker-ce/ubuntu/) for details.
3. Shut down the VM, expand the VM disk to desired size and restart the VM.
4. download the reference data (reference genome, wig files etc) to `/refdata` directory on the VM.
5. run ```sudo waagent -deprovision+user``` on the VM
6. go to the azure portal and capture the vm. 
7. the ID can be found under the ```Resource ID``` field of the image.

NOTE: batch only supports ubuntu 16.04 at the moment, images built with 17.04 will not work.

