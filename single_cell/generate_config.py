'''
Created on Apr 9, 2018

@author: dgrewal
'''

from single_cell.config import pipeline_config
from single_cell.config import batch


def generate_config(args):
    config_yaml = args.get("pipeline_config")
    batch_yaml = args.get("batch_config")
    params_override = args.get("config_override")

    if config_yaml:
        config_params = pipeline_config.get_config_params(override=params_override)
        config = pipeline_config.get_singlecell_pipeline_config(config_params)
        pipeline_config.write_config(config, config_yaml)

    if batch_yaml:
        config_params = batch.get_batch_params(override=params_override)
        config = batch.get_batch_config(config_params)
        batch.write_config(config, batch_yaml)
