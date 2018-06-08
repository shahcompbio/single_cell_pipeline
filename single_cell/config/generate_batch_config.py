
import os
from single_cell.utils import helpers
import warnings
import batch


def generate_submit_config_in_temp(args):

    if args['which'] in ['clean_sentinels', 'generate_config']:
        return args

    if args.get("submit_config", None):
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
        warnings.warn("no tmpdir specified, generating configs in working dir")
        batch_yaml = os.path.join(os.getcwd(), batch_yaml)

    batch_yaml = helpers.get_incrementing_filename(batch_yaml)

    params_override = args["config_override"]

    config_params = batch.get_batch_params(override=params_override)
    config = batch.get_batch_config(config_params)
    batch.write_config(config, batch_yaml)

    args["submit_config"] = batch_yaml

    return args
