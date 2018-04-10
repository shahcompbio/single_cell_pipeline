'''
Created on Apr 9, 2018

@author: dgrewal
'''


from cmdline import parse_args

from run import run_pipeline
from generate_config import generate_config

def main():

    args = parse_args()

    if "generate_config" in args:
        generate_config(args)

    else:
        run_pipeline(args)




if __name__ == "__main__":
    main()