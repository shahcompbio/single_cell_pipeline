import os
import logging
from single_cell.config import pipeline_config
from single_cell.utils import helpers


def generate_pipeline_config_in_temp(args):

    if args['which'] in ['clean_sentinels', 'generate_config']:
        return args

    if args.get("config_file", None):
        return args

    config_yaml = "config.yaml"
    tmpdir = args.get("tmpdir", None)
    pipelinedir = args.get("pipelinedir", None)

    # use pypeliner tmpdir to store yaml
    if pipelinedir:
        config_yaml = os.path.join(pipelinedir, config_yaml)
    elif tmpdir:
        config_yaml = os.path.join(tmpdir, config_yaml)
    else:
        logging.getLogger("single_cell.generate_pipeline_config").warn(
            "no tmpdir specified, generating configs in working dir"
        )
        config_yaml = os.path.join(os.getcwd(), config_yaml)

    config_yaml = helpers.get_incrementing_filename(config_yaml)

    params_override = args["config_override"]

    helpers.makedirs(config_yaml, isfile=True)

    config_params = pipeline_config.get_config_params(override=params_override)
    config = pipeline_config.get_singlecell_pipeline_config(config_params, override=params_override)
    pipeline_config.write_config(config, config_yaml)

    args["config_file"] = config_yaml

    return args
