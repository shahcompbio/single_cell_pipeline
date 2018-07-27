# Single Cell Pipeline on Azure

## Azure Blob
1. Install all dependencies. See ```build_vm.sh``` for details. 
2. Setup an Azure Blob Storage account.
3. setup Azure Active Directory. Azure AD should have write permissions over storage account. see https://docs.microsoft.com/en-us/azure/azure-resource-manager/resource-group-create-service-principal-portal for details.
4. source:
	```
	export CLIENT_ID='<Active-Directory-client-id>'
	export TENANT_ID='<Active-Directory-tenant-id>'
	export SECRET_KEY='<Active-Directory-secret>'
	export SUBSCRIPTION_ID='<Azure-subscription-id>'
	export RESOURCE_GROUP='<resource group name>'
	```
5. Upload data to blob storage
6. build the input file. The overall structure of the file stays the same except for file paths. The file paths follow ```<storage_account_name>/<container_name>/<blob path>``` structure. Same applies for the ```--pipelinedir``` , ```--tmpdir```and ```--output``` directories
7. run the pipeline. add ```--storage azureblob``` to the command to use blob backend.

NOTE: The pipeline supports multiple Azure storage accounts. please ensure all storage accounts are in the same resource group and that azure AD has write permissions over all groups. To use multiple accounts, simply replace storage_account_name in your file path with the relevant account name.


## Azure Batch
The Azure Batch backend requires the Azure Blob, local storage or any other kind of storage is not supported in this mode.

1. Start a small vm. this will serve as our master node. Install singlecellpipeline and pypeliner.
2. create a new batch account. see https://docs.microsoft.com/en-us/azure/batch/batch-account-create-portal for details
3. Create a new storage account (or use an existing one).
4. setup Azure Active Directory. The batch account and storage account(s) must be in the same resource group.  Azure AD should have write permissions over batch account and the storage account(s). see https://docs.microsoft.com/en-us/azure/azure-resource-manager/resource-group-create-service-principal-portal for details. 
5. Set up the image. see ```Build Azure Managed Image``` for details.
5. source on headnode:
	```
	export AZURE_BATCH_URL='<batch-url>'
	export CLIENT_ID='<Active-Directory-client-id>'
	export TENANT_ID='<Active-Directory-tenant-id>'
	export SECRET_KEY='<Active-Directory-secret>'
	export SUBSCRIPTION_ID='<Azure-subscription-id>'
	export RESOURCE_GROUP='<resource group name>'
	export IMAGE_ID='<image-to-start-vms-from>'
	export SKU='batch.node.ubuntu 16.04'
	``` 
6.  Launch the pipeline. add ``` --submit azurebatch``` to use the batch backend.

## Rabbit MQ
The Azure Storage Accounts have a limited amount of available bandwidth. see  https://docs.microsoft.com/en-us/azure/storage/common/storage-scalability-targets for details. 
Unfortunately Azure Storage Accounts do not  throttle the bandwidth when the storage accounts go over the scalability limits, instead all open connections die with an ```AzureHttpError```.  The single cell pipeline supports RabbitMq messaging queue to throttle the number of open connections to avoid these errors.
This part requires a server running RabbitMQ. see https://docs.bitnami.com/azure/infrastructure/rabbitmq/ to start a rabbitmq server on azure.

Please export the following environment variables to enable rabbitmq throttling:
```
export RABBITMQ_USERNAME="<username>"
export RABBITMQ_PASSWORD="<password>"
export RABBITMQ_IP="<ip address>"
export RABBITMQ_VHOST="<vhostname>"
```
Please ensure the vhost is available prior to launching the pipeline. The pipeline will start one queue per storage account under the vhost. 


## Build Azure Managed Image

1. start a ubuntu 16.04 vm. 
2. Install all dependencies. See ```build_vm.sh``` for details.
3. attach 2 managed disks. A smaller (under 64gb) disk for storing reference data and a bigger disk (preferably 1 TB). mount small disk at /refdata and the bigger one at /datadrive.
4. download the reference data (reference genome, wig files etc) to the reference disk.
5. run ```sudo waagent -deprovision+user```
6. go to the azure portal and capture the vm. 
7. the image id can be found under the ```Resource ID``` field

NOTE: batch only supports ubuntu 16.04 at the moment, images built with 17.04 will not work.

