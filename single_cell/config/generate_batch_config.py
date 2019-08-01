
import logging
import os

from single_cell.config import batch
from single_cell.utils import helpers


def generate_submit_config_in_temp(args):

    if args['which'] in ['clean_sentinels', 'generate_config']:
        return args

    if args.get("submit_config", None):
        return args

    azure_submit = ['azurebatch',
                    'pypeliner.contrib.azure.batchqueue.AzureJobQueue']
    if not args.get("submit", None) in azure_submit:
        return args

    batch_yaml = "batch.yaml"
    tmpdir = args.get("tmpdir", None)
    pipelinedir = args.get("pipelinedir", None)

    # use pypeliner tmpdir to store yaml
    if pipelinedir:
        batch_yaml = os.path.join(pipelinedir, batch_yaml)
    elif tmpdir:
        batch_yaml = os.path.join(tmpdir, batch_yaml)
    else:
        logging.getLogger("single_cell.generate_batch_config").warn(
            "no tmpdir specified, generating configs in working dir"
        )
        batch_yaml = os.path.join(os.getcwd(), batch_yaml)

    helpers.makedirs(batch_yaml, isfile=True)

    batch_yaml = helpers.get_incrementing_filename(batch_yaml)

    params_override = args["config_override"]
    if not params_override:
        params_override = {}

    config_params = batch.get_batch_params(override=params_override)
    config = batch.get_batch_config(config_params, override=params_override)
    batch.write_config(config, batch_yaml)

    args["submit_config"] = batch_yaml

    return args
