import argparse
import os

import yaml
from single_cell.utils import helpers
from single_cell.utils import storageutils

import compare


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


def get_storage_account(path):
    path = path.strip()
    return path.split('/')[0]


def download_blob(blob_path, tempdir):
    outpath = os.path.join(tempdir, blob_path)
    helpers.makedirs(outpath, isfile=True)

    storageutils.download_blob(blob_path, outpath, storage='azureblob')

    return outpath


def compare_output(reads, metrics, ref_reads, ref_metrics, tempdir):
    reads = download_blob(reads, tempdir)
    metrics = download_blob(metrics, tempdir)
    ref_reads = download_blob(ref_reads, tempdir)
    ref_metrics = download_blob(ref_metrics, tempdir)

    compare.compare_metrics(ref_metrics, metrics)
    compare.compare_reads(ref_reads, reads)


def parse_args():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers()

    container_yaml = subparsers.add_parser("generate_container_yaml")
    container_yaml.set_defaults(which='container_yaml')
    container_yaml.add_argument('output',
                                help='specify path to the output file')

    compare = subparsers.add_parser('compare')
    compare.set_defaults(which='compare')
    compare.add_argument(
        '--ref_reads',
        help='specify path to the output dir'
    )
    compare.add_argument(
        '--ref_metrics',
        help='specify path to the output dir'
    )
    compare.add_argument(
        '--reads',
        help='specify path to the output dir'
    )
    compare.add_argument(
        '--metrics',
        help='specify path to the output dir'
    )
    compare.add_argument(
        '--tempdir',
        help='specify path to the output dir'
    )

    args = parser.parse_args()

    args = vars(args)
    return args


def main(args):
    if args['which'] == 'container_yaml':
        generate_container_yaml(args['output'])

    elif args['which'] == 'compare':
        compare_output(
            args['reads'], args['metrics'],
            args['ref_reads'], args['ref_metrics'],
            args['tempdir']
        )


if __name__ == '__main__':
    args = parse_args()
    main(args)
