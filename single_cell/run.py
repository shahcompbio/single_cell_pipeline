import sys

from single_cell.aneufinder import aneufinder_pipeline
from single_cell.breakpoint_calling import breakpoint_calling_pipeline
from single_cell.clean_sentinels import clean_sentinels
from single_cell.cmdline import parse_args
from single_cell.copy_number import copy_number_calling_pipeline
from single_cell.docker_run import run_with_docker
from single_cell.generate_config import generate_config
from single_cell.germline_calling import germline_calling_pipeline
from single_cell.infer_haps import infer_haps_pipeline
from single_cell.ltm import ltm_pipeline
from single_cell.merge_bams import merge_bams_pipeline
from single_cell.multi_sample import multi_sample_pipeline
from single_cell.qc import qc_pipeline
from single_cell.split_bam import split_bam_pipeline
from single_cell.variant_calling import variant_calling_pipeline
from single_cell.variant_calling import variant_counting_pipeline


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
