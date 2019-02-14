# Pypeliner on AWS

### Build an AMI
1. Go to Services -> EC2 -> Launch Instance
2. *Choose AMI*: Pick the latest Amazon Linux 2 AMI at the AMI page
3. *Instance Size*: Pick a VM size, t2.micro should be enough
4. Choose *Review and Launch*.
5. *Review Instance Launch*: Go to storage and edit. Change the size of the disk from 8GB to 512GB. Change from SSD to Magnetic.
6. Specify your Key-Value pair if you're asked for it, create a new key value pair if you dont already have one.
7. go to Services -> EC2 -> Running Instances. Track down your VM and copy the IPv4 public IP address. 
8. Log into the VM with `ssh ec2-user@<ip-address-here>`
9. Install ecs-init on the node. `sudo amazon-linux-extras install -y ecs; sudo systemctl enable --now ecs`. If you get `amazon-linux-extras: command not found` error, then your EC2 instance is not running Amazon Linux 2 AMI, chances are you're running Amazon Linux AMI.
10. create the refdata directory and upload the reference files. `mkdir /refdata` and then rsync/download data from these locations: 
    a) Download from Azure. The data is under singlecelldata storage account and in reference container.
    b) Ask Diljot <grewald@mskcc.org> for access to S3 bucket with the data.
11. Go to Services -> EC2 -> Running Instances and track down your VM. Right click to open the menu and choose Image -> Create Image. Specify the image name and choose Create. The image should show up in EC2 -> AMIs. locate the AMI-id under the details pane, you'll need this later.

Detailed instructions for starting a new EC2 instance at: [Starting AWS EC2 instance](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/EC2_GetStarted.html#ec2-launch-instance).

### Start Master EC2 instance.
1. Go to EC2 -> AMIs. Pick the AMI we created in previous section. Choose Launch from the drop down menu.
2. Choose an instance size. t2.large would be a good starting point. Choose Review and Launch
3. Choose Launch.
4. Specify your Key-Value pair if you're asked for it, create a new key value pair if you dont already have one.
5. go to Services -> EC2 -> Running Instances. Track down your VM and copy the IPv4 public IP address. 
6. Log into the VM with `ssh ec2-user@<ip-address-here>`
7. You need access to docker containers required to run the pipeline next. The dockerfiles are at [github](https://github.com/shahcompbio/docker_containers). You can either
    a) build the containers and push to ECR
    b) Ask Diljot <grewald@mskcc.org> for access to ECR containers
8. Install AWS CLI. see [AWS](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-install.html)
9. Run `aws configure`
10. log into AWS ECR by running `aws ecr get-login --no-include-email`
11. go to Services -> ECR and get your registry server name. The name likely follows this format: `<number>.dkr.ecr.us-east-2.amazonaws.com`.
12. pull the latest single_cell_pipeline docker image with `docker pull <registry-server-name>/scp/single_cell_pipeline:<version>`
13. You should now have a working environment on the machine and can run the pipeline with: `docker run <registry-server-name>/scp/single_cell_pipeline:<version> single_cell -h`


### Run on AWS batch


