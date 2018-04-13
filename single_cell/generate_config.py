'''
Created on Apr 9, 2018

@author: dgrewal
'''

import os
from single_cell.config import generate_pipeline_config
from single_cell.config import generate_batch_config


def generate_config(args):
    config_yaml = args["pipeline_config"]
    batch_yaml = args["batch_config"]
    params_override = args["config_override"]

    generate_pipeline_config.main(
        output=config_yaml,
        input_params=params_override)

    generate_batch_config.main(output=batch_yaml, input_params=params_override)
