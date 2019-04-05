# Set up Azure Batch


### 1. Create a resource group

### 2. Create a new storage account
Please follow the steps at [azure-storage-quickstart](https://docs.microsoft.com/en-us/azure/storage/common/storage-quickstart-create-account?tabs=azure-portal)

### 3. Create a new batch account
Please follow the steps at [azure-batch-quickstart](https://docs.microsoft.com/en-us/azure/batch/batch-account-create-portal)

### 4. Create a new docker container registry (optional)
Please follow the steps at [azure-docker-registry-quickstart](https://docs.microsoft.com/en-us/azure/container-registry/container-registry-get-started-portal)

### 5. Create a new keyvault
Please follow the steps at [azure-keyvault-quickstart](https://docs.microsoft.com/en-us/azure/key-vault/quick-create-portal)

### 6. Create a new Active directory app.

1. Please follow the steps at [create-service-principal-for-azure-ad](https://docs.microsoft.com/en-us/azure/azure-stack/user/azure-stack-create-service-principals#create-service-principal-for-azure-ad). After this step, you should have an application ID for the service principal and a secret key.
2. Please follow the steps at [use-an-existing-tenant](https://docs.microsoft.com/en-us/azure/active-directory/develop/quickstart-create-new-tenant#use-an-existing-tenant) to retrieve your tenant ID.
3. you can also find the subscription ID under the subscriptions page from the portal.

### 7. set the service principal role:
1. Go to the resource group we created in step 1.
2. From the left pane, select *Access Control(IAM)*.
3. Go to role assignments and add a new role assignment. Choose the AD service principal from step 6 and give it owner permissions.

### 8. set up keyvault for the storage account
1. go to the storage account from step 2. choose Access Keys from the left pane. Copy the key1 value.
2. go to the keyvault from step 5. choose secrets and then create a new secret. Set the storage account name as the name and paste the storage account key in value.
3. Choose Access Policies from the left pane of the keyvault. Choose Add New, select Secret Management profile and choose service principal under the authorized application. click Ok and then save.

### 9. Create an image:
Now we need to build a new Azure Image. This image will be used start the nodes that will run the jobs in batch. So this image should be built from a node that contains all required dependencies.

#### Create a new VM
Please follow [create-vm](https://docs.microsoft.com/en-us/azure/virtual-machines/linux/quick-create-portal) to start a new linux VM in the resource group from 1. Please make sure that its an ubuntu 18.04 machine.

#### Install pipeline dependencies
###### 1. docker
If the pipeline that you intend to run is configured to run docker then you can just install docker on the image. Please follow [install-docker](https://docs.docker.com/install/linux/docker-ce/ubuntu/#install-docker-ce-1) to install docker on the node
###### 2. no docker
If you're not using docker then please go ahead and install all dependencies on the node with/without the package manager of your choice. The VM should have exactly the same environment as the node that will be used to launch the workflow.

NOTE: Every change in the dependencies will require a new image which can be very time consuming. We recommend docker for precisely this reason

#### Download reference data
Download the reference data to a location of your choice on the node. 

NOTE: Please feel free to add more disks or resize the main OS disk. We recommend resizing the main OS disk to make the setup easier. Please follow [expand-disks](https://docs.microsoft.com/en-us/azure/virtual-machines/linux/expand-disks)

Follow [capture-image-resource](https://docs.microsoft.com/en-us/azure/virtual-machines/windows/capture-image-resource) to capture an image from this node.
