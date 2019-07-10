import argparse
import os
import errno

import azure.storage.blob as azureblob
import yaml
from azure.common.credentials import ServicePrincipalCredentials
from azure.keyvault import KeyVaultClient, KeyVaultAuthentication
from azure.keyvault.models import KeyVaultErrorException


def makedirs(directory, isfile=False):

    if isfile:
        directory = os.path.dirname(directory)
        if not directory:
            return

    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


class UnconfiguredStorageAccountError(Exception):
    pass


def generate_container_yaml(filepath):
    data = {
        'docker': {
            'server': os.environ['CONTAINER_REGISTRY'],
            'username': os.environ['CLIENT_ID'],
            'password': os.environ['SECRET_KEY'],
            'mounts': {
                'refdata': '/refdata',
                'mnt': '/mnt',
                'datadrive': '/datadrive'
            }
        }
    }

    with open(filepath, 'w') as container_out:
        yaml.dump(data, container_out)


def get_storage_account_key(
        accountname, client_id, secret_key, tenant_id, keyvault_account
):
    """
    Uses the azure management package and the active directory
    credentials to fetch the authentication key for a storage
    account from azure key vault. The key must be stored in azure
    keyvault for this to work.
    :param str accountname: storage account name
    """

    def auth_callback(server, resource, scope):
        credentials = ServicePrincipalCredentials(
            client_id=client_id,
            secret=secret_key,
            tenant=tenant_id,
            resource="https://vault.azure.net"
        )
        token = credentials.token
        return token['token_type'], token['access_token']

    client = KeyVaultClient(KeyVaultAuthentication(auth_callback))
    keyvault = "https://{}.vault.azure.net/".format(keyvault_account)
    # passing in empty string for version returns latest key
    try:
        secret_bundle = client.get_secret(keyvault, accountname, "")
    except KeyVaultErrorException:
        err_str = "The pipeline is not setup to use the {} account. ".format(accountname)
        err_str += "please add the storage key for the account to {} ".format(keyvault_account)
        err_str += "as a secret. All input/output paths should start with accountname"
        raise UnconfiguredStorageAccountError(err_str)
    account_key = secret_bundle.value

    return account_key


def get_blob_client(storage_account_name):
    storage_account_key = get_storage_account_key(
        storage_account_name,
        os.environ['CLIENT_ID'],
        os.environ['SECRET_KEY'],
        os.environ['TENANT_ID'],
        os.environ['AZURE_KEYVAULT_ACCOUNT'],
    )

    blob_client = azureblob.BlockBlobService(
        account_name=storage_account_name,
        account_key=storage_account_key)
    return blob_client


def download_blob(blob_client, container, blobpath, filepath):
    blob = blob_client.get_blob_to_path(
        container,
        blobpath,
        filepath)
    return blob


def list_all_blobs(blob_client, container):
    blobs = blob_client.list_blobs(container)
    for blob in blobs:
        yield blob.name


def download_refdata(basedir):
    blob_client = get_blob_client(os.environ['REFDATA_STORAGE_ACCOUNT'])
    for blob in list_all_blobs(blob_client, 'refdata'):
        outpath = os.path.join(basedir, blob)
        if not os.path.exists(outpath):
            makedirs(outpath, isfile=True)
            print "downloading: {}".format(blob)
            download_blob(blob_client, 'refdata', blob, outpath)


def run_pipeline():
    pass


def parse_args():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers()

    container_yaml = subparsers.add_parser("generate_container_yaml")
    container_yaml.set_defaults(which='container_yaml')
    container_yaml.add_argument('output',
                                help='specify path to the output file')

    download_refdata = subparsers.add_parser('download_refdata')
    download_refdata.set_defaults(which='download_refdata')
    download_refdata.add_argument('basedir',
                                help='specify path to the output dir')

    args = parser.parse_args()

    args = vars(args)
    return args


def main(args):
    if args['which'] == 'container_yaml':
        generate_container_yaml(args['output'])
    elif args['which'] == 'download_refdata':
        download_refdata(args['basedir'])


if __name__ == '__main__':
    args = parse_args()
    main(args)

