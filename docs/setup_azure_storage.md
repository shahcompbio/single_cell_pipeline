
## Blob Storage
In this section, we will setup single cell pipeline to use Azure Blob Storage. Each job will download the data from a blob storage account and push the outputs to blob storage.

### 1. Create a Storage account.
Create/Enter a resource group -> Add -> Storage account - blob, file, table, queue. \
See [link](https://docs.microsoft.com/en-us/azure/storage/common/storage-account-create?tabs=azure-portal) for documentation
### 2. Upload data to Blob Storage.
  1. Use portal:
  Enter a storage account -> Tools and SDKs (Storage Explorer) -> Go to the destination container -> Upload/Delete/New Folder/etc..

  2. Use Azure Cli:

  ```
  az storage blob upload-batch -d container_name/prefix -s path/to/data/to/upload --acount-name destination_storage_account_name --pattern some_pattern_here --auth_mode login
  ```
*Reference:*\
*Manage storage account and data using Azure Portal:*  [Azure Portal Quickstart guide](https://docs.microsoft.com/en-us/azure/storage/blobs/storage-quickstart-blobs-portal)\
*Manage a storage account using Azure Cli:*  [Official Documentation](https://docs.microsoft.com/en-us/azure/storage/blobs/storage-quickstart-blobs-cli)\
*Manage data in a storage account using Azure Cli:* [az storage blob documentation](https://docs.microsoft.com/en-us/cli/azure/storage/blob?view=azure-cli-latest#az-storage-blob-upload-batch)

### 3. Create another Storage account for pypeliner temp files.
(You can also use the same storage account as step 1)\
You would need to create following temporary directories to store the temp files of the pipeline.
```
temp
output
bam
pipeline
```
### 4. Setup an azure key vault.

1. Add a Key Vault to your Resource Group. Please follow "Create a vault" section at [Link](https://docs.microsoft.com/en-us/azure/key-vault/quick-create-portal)
2. Add a secret to the KeyVault. Please follow "Add a secret to Key Vault" section at [Link](https://docs.microsoft.com/en-us/azure/key-vault/quick-create-portal). Please use the storage account name as the key and the storage account key as the value. The storage account key can be accessed from storage account -> access keys -> key1
3. repeat step 2 for all storage accounts used by the pipeline.

### 5. Setup an active directory app.

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


### 6. Launch the pipeline

1. Please add `--storage azureblob` to use the azure blob storage backend.

  Example command using conda environment and Azure Blob storage:
  ```
  single_cell alignment \
  --input_yaml /path/to/inputs.yaml \
  --library_id A97318A \
  --maxjobs 4 \
  --nocleanup \
  --sentinel_only  \
  --config_override '{"refdir": "storage_account/container/refdata"}' \
  --submit local \
  --loglevel DEBUG \
  --tmpdir storage_account/container/temp \
  --pipelinedir storage_account/container/pipeline \
  --submit local \
  --out_dir storage_account/container/output \
  --bams_dir storage_account/container/bams

  ```
Example command using docker and Azure Blob storage (Make sure you put the environmental vars for Azure in a file and provide the path to it):
```
docker run -w /datadrive -v /home:/home -v /datadrive:/datadrive -v /var/run/docker.sock:/var/run/docker.sock -v /usr/bin/docker:/usr/bin/docker
 --env-file /path/to/file/contains/env/vars \
 --rm singlecellpipeline/single_cell_pipeline:v0.5.17 single_cell alignment \
 --input_yaml /path/to/inputs.yaml \
 --storage azureblob \
 --library_id A97318A \
 --maxjobs 1 \
 --nocleanup \
 --sentinel_only  \
 --config_override '{"refdir": "/path/to/refdata"}' \
 --context_config /path/to/context_config.yaml \
 --loglevel DEBUG \
 --tmpdir storage_account/container/temp \
 --pipelinedir storage_account/container/pipeline \
 --submit local \
 --out_dir storage_account/container/output \
 --bams_dir storage_account/container/bams
```
2. All input files should be uploaded to blob prior to launching the pipeline.
3. The file paths should be formatted as `{storage-account-name}/{container-name}/{blob-path}`. Paths to the `--pipelinedir` and `--tmpdir` should also follow the same format.

  Example ```inputs.yaml```:
  ```
SA1090-A96213A-R20-C28:
  column: 28
  condition: B
  fastqs:
    HHCJ7CCXY_5.HGTJJCCXY_8.HYG5LCCXY_6.HYG5LCCXY_7.HYG5LCCXY_5:
      fastq_1: scpipelinetestdata/data/SA1090-A96213A-R20-C28_1.fastq.gz
      fastq_2: scpipelinetestdata/data/SA1090-A96213A-R20-C28_2.fastq.gz
      sequencing_center: TEST
      sequencing_instrument: TEST
  img_col: 45
  index_i5: i5-20
  index_i7: i7-28
  pick_met: C1
  primer_i5: GTATAG
  primer_i7: CTATCT
  row: 20
SA1090-A96213A-R20-C62:
  column: 62
  condition: B
  fastqs:
    HHCJ7CCXY_5.HGTJJCCXY_8.HYG5LCCXY_6.HYG5LCCXY_7.HYG5LCCXY_5:
      fastq_1: scpipelinetestdata/data/SA1090-A96213A-R20-C62_1.fastq.gz
      fastq_2: scpipelinetestdata/data/SA1090-A96213A-R20-C62_2.fastq.gz
      sequencing_center: TEST
      sequencing_instrument: TEST
  img_col: 11
  index_i5: i5-20
  index_i7: i7-62
  pick_met: C1
  primer_i5: GTATAG
  primer_i7: AAGCTA
  row: 20
SA1090-A96213A-R22-C43:
  column: 43
  condition: B
  fastqs:
    HHCJ7CCXY_5.HGTJJCCXY_8.HYG5LCCXY_6.HYG5LCCXY_7.HYG5LCCXY_5:
      fastq_1: scpipelinetestdata/data/SA1090-A96213A-R22-C43_1.fastq.gz
      fastq_2: scpipelinetestdata/data/SA1090-A96213A-R22-C43_2.fastq.gz
      sequencing_center: TEST
      sequencing_instrument: TEST
  img_col: 30
  index_i5: i5-22
  index_i7: i7-43
  pick_met: C2
  primer_i5: GCTGTA
  primer_i7: ATTCCG
  row: 22
  ```
