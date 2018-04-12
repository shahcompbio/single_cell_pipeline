'''
Created on Apr 9, 2018

@author: dgrewal
'''

import os
from single_cell.config import config_generator


def generate_config(args):


    cfgdir = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'config')

    output = args["output"]

    input_params = args["config_override"]

    config = os.path.join(cfgdir, "config.yaml")

    params = os.path.join(cfgdir, "params.yaml")


    config_generator.main(output=output, input_params = input_params)
