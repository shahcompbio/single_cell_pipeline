import sys

from aneufinder import aneufinder_pipeline
from breakpoint_calling import breakpoint_calling_pipeline
from clean_sentinels import clean_sentinels
from cmdline import parse_args
from copy_number import copy_number_calling_pipeline
from docker_run import run_with_docker
from generate_config import generate_config
from germline_calling import germline_calling_pipeline
from infer_haps import infer_haps_pipeline
from ltm import ltm_pipeline
from merge_bams import merge_bams_pipeline
from multi_sample import multi_sample_pipeline
from qc import qc_pipeline
from split_bam import split_bam_pipeline
from variant_calling import variant_calling_pipeline
from variant_calling import variant_counting_pipeline


def main():
    args = parse_args()

    if args["which"] == "generate_config":
        generate_config(args)
        return

    if args["which"] == "clean_sentinels":
        clean_sentinels(args)
        return

    if args["run_with_docker"]:
        run_with_docker(args, sys.argv)
        return

    if args["which"] == "qc":
        qc_pipeline(args)

    if args["which"] == "aneufinder":
        aneufinder_pipeline(args)

    if args["which"] == "merge_bams":
        merge_bams_pipeline(args)

    if args["which"] == "split_bam":
        split_bam_pipeline(args)

    if args["which"] == "variant_calling":
        variant_calling_pipeline(args)

    if args["which"] == "copy_number_calling":
        copy_number_calling_pipeline(args)

    if args["which"] == "infer_haps":
        infer_haps_pipeline(args)

    if args["which"] == "germline_calling":
        germline_calling_pipeline(args)

    if args["which"] == "breakpoint_calling":
        breakpoint_calling_pipeline(args)

    if args["which"] == "variant_counting":
        variant_counting_pipeline(args)

    if args["which"] == "ltm":
        ltm_pipeline(args)

    if args["which"] == "multi_sample_pseudo_bulk":
        multi_sample_pipeline(args)


if __name__ == "__main__":
    main()
