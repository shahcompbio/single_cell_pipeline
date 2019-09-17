import os
import shutil
from pypeliner.contrib.azure import blobclient


class UnknownStorage(Exception):
    pass


def get_storage_account(path):
    path = path.strip()
    return path.split('/')[0]


def download_azure_blob(blob_path, local_path):
    client_id = os.environ["CLIENT_ID"]
    secret_key = os.environ["SECRET_KEY"]
    tenant_id = os.environ["TENANT_ID"]
    keyvault_account = os.environ['AZURE_KEYVAULT_ACCOUNT']

    try:
        os.makedirs(os.path.dirname(local_path))
    except OSError:
        pass

    storageaccountname = get_storage_account(blob_path)

    client = blobclient.BlobStorageClient(
        storageaccountname, client_id, tenant_id, secret_key,
        keyvault_account
    )

    client.download_to_path(local_path, blob_uri=blob_path)


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
