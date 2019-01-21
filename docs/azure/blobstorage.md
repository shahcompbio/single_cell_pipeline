
## Blob Storage
In this section, we will setup pypeliner to use Azure Blob Storage. Each job will download the data from a blob storage account and push the outputs to blob storage. 

1. Create a Storage account. see [link](https://code.visualstudio.com/tutorials/static-website/create-storage) for documentation
2. Upload data to Blob Storage. see [portal](https://docs.microsoft.com/en-us/azure/storage/blobs/storage-quickstart-blobs-portal) and [CLI](https://docs.microsoft.com/en-us/azure/storage/blobs/storage-quickstart-blobs-cli)
3. Create another Storage account for pypeliner temp files. You can also use the same storage account as step 1
4. setup an azure key vault.
5. setup an active directory app. 
6. Run the pipeline 

### Key Vault
1. Add a Key Vault to your Account. Please follow "Create a vault" section at [Link](https://docs.microsoft.com/en-us/azure/key-vault/quick-create-portal)
2. Add a secret to the KeyVault. Please follow "Add a secret to Key Vault" section at [Link](https://docs.microsoft.com/en-us/azure/key-vault/quick-create-portal). Please use the storage account name as the key and the storage account key as the value. The storage account key can be accessed from storage account -> access keys -> key1
3. repeat step 2 for all storage accounts used by the pipeline.
### Active Directory
1. set up an Active directory app. see [link](https://docs.microsoft.com/en-us/azure/active-directory/develop/app-objects-and-service-principals) for documentation.
2. source the following env variables
    ```
    export AZURE_BATCH_URL='https://scbatch.canadaeast.batch.azure.com'
    export AZURE_STORAGE_ACCOUNT=''
    export CLIENT_ID=''
    export TENANT_ID=''
    export SECRET_KEY=''
    export SUBSCRIPTION_ID=''
    export RESOURCE_GROUP=''
    export AZURE_KEYVAULT_ACCOUNT=''
    ```
3. Go to your resource group -> Access Control -> Add Role Assignment. Give 'owner' access to Active Directory app from step 1.
### Launch the pipeline
1. Please add `--storage azureblob` to use the azure blob storage backend.
2. All input files should be uploaded to blob prior to launching the pipeline. 
3. The paths should be formatted as `/{storage-account-name}/{container-name}/{blob-path}`. Paths to the `--pipelinedir` and `--tmpdir` should also follow the same format
 
### Rabbit MQ
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

