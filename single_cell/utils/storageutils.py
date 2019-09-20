import os
import shutil
from pypeliner.contrib.azure import blobclient
from single_cell.utils  import helpers

class UnknownStorage(Exception):
    pass


def get_storage_account(path):
    path = path.strip('/')
    return path.split('/')[0]

def unpack_path(path):
    path = path.strip('/').split('/')

    storage_account = path[0]
    container = path[1]
    blob_name = '/'.join(path[2:])

    return storage_account, container, blob_name


def download_azure_blob(blob_path, local_path):
    helpers.makedirs(local_path, isfile=True)

    client_id = os.environ["CLIENT_ID"]
    secret_key = os.environ["SECRET_KEY"]
    tenant_id = os.environ["TENANT_ID"]
    keyvault_account = os.environ['AZURE_KEYVAULT_ACCOUNT']

    storageaccountname = get_storage_account(blob_path)

    client = blobclient.BlobStorageClient(
        storageaccountname, client_id, tenant_id, secret_key,
        keyvault_account
    )

    client.download_to_path(local_path, blob_uri=blob_path)


def download_azure_blobs(prefix, local_dir):
    helpers.makedirs(local_dir)

    client_id = os.environ["CLIENT_ID"]
    secret_key = os.environ["SECRET_KEY"]
    tenant_id = os.environ["TENANT_ID"]
    keyvault_account = os.environ['AZURE_KEYVAULT_ACCOUNT']

    storage_account, container, blob_prefix = unpack_path(prefix)

    client = blobclient.BlobStorageClient(
        storage_account, client_id, tenant_id, secret_key,
        keyvault_account
    )

    for blob_uri in client.list_blobs(container_name=container, prefix=blob_prefix):
        download_to_path = os.path.join(local_dir, blob_uri.name)
        helpers.makedirs(download_to_path, isfile=True)

        download_uri = os.path.join(storage_account, container, blob_uri.name)

        client.download_to_path(download_to_path, blob_uri=download_uri)



def upload_azure_blob(blob_path, filepath):
    client_id = os.environ["CLIENT_ID"]
    secret_key = os.environ["SECRET_KEY"]
    tenant_id = os.environ["TENANT_ID"]
    keyvault_account = os.environ['AZURE_KEYVAULT_ACCOUNT']

    storageaccountname = get_storage_account(blob_path)

    client = blobclient.BlobStorageClient(
        storageaccountname, client_id, tenant_id, secret_key,
        keyvault_account
    )

    client.upload_from_file(
        filepath, blob_uri=blob_path
    )


def upload_azure_blobs(blob_prefix, local_dir):
    raise NotImplementedError()


def download_aws_blob(blob_path, local_path):
    raise NotImplementedError()


def upload_aws_blob(blob_path, local_path):
    raise NotImplementedError()


def download_blob(blob_path, local_path, storage=None):
    if not storage:
        shutil.move(blob_path, local_path)
    elif storage == 'azureblob':
        download_azure_blob(blob_path, local_path)
    elif storage == 'awss3':
        download_aws_blob(blob_path, local_path)
    else:
        raise UnknownStorage(storage)


def upload_blob(blob_path, local_path, storage=None):
    if not storage:
        shutil.move(local_path, blob_path)
    elif storage == 'azureblob':
        upload_azure_blob(blob_path, local_path)
    elif storage == 'awss3':
        upload_aws_blob(blob_path, local_path)
    else:
        raise UnknownStorage(storage)


def download_blobs(blob_prefix, local_dir, storage=None):
    if not storage:
        shutil.move(blob_prefix, local_dir)
    elif storage == 'azureblob':
        download_azure_blobs(blob_prefix, local_dir)
    elif storage == 'awss3':
        download_aws_blob(blob_prefix, local_dir)
    else:
        raise UnknownStorage(storage)
