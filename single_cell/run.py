import pypeliner
from cmdline import parse_args
from align import align_workflow
from hmmcopy import hmmcopy_workflow
from aneufinder import aneufinder_workflow
from merge_bams import merge_bams_workflow
from variant_calling import variant_calling_workflow
from variant_calling import variant_counting_workflow
from germline_calling import germline_calling_workflow
from breakpoint_calling import breakpoint_calling_workflow
from split_bam import split_bam_workflow
from generate_config import generate_config
from clean_sentinels import clean_sentinels
from copy_number import copy_number_calling_workflow
from ltm import ltm_workflow
from infer_haps import infer_haps_workflow
from multi_sample import multi_sample_workflow
# from copyclone import copyclone_workflow
from docker_run import run_with_docker
import sys


def main():

    args = parse_args()

    if args["run_with_docker"]:
        run_with_docker(args, sys.argv)
        return

    if args["which"] == "generate_config":
        generate_config(args)
        return

    if args["which"] == "clean_sentinels":
        clean_sentinels(args)
        return

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = None

    if args["which"] == "align":
        workflow = align_workflow(args)

    if args["which"] == "hmmcopy":
        workflow = hmmcopy_workflow(args)

#     if args["which"] == "copyclone":
#         pyp = pypeliner.app.Pypeline(config=args)
#         workflow = pypeliner.workflow.Workflow()
#         workflow = copyclone_workflow(workflow, args)
#         pyp.run(workflow)

    if args["which"] == "aneufinder":
        workflow = aneufinder_workflow(args)

    if args["which"] == "merge_bams":
        workflow = merge_bams_workflow(args)

    if args["which"] == "split_bam":
        workflow = split_bam_workflow(args)

    if args["which"] == "variant_calling":
        workflow = variant_calling_workflow(args)

    if args["which"] == "copy_number_calling":
        workflow = copy_number_calling_workflow(args)

    if args["which"] == "infer_haps":
        workflow = infer_haps_workflow(args)

    if args["which"] == "germline_calling":
        workflow = germline_calling_workflow(args)

    if args["which"] == "breakpoint_calling":
        workflow = breakpoint_calling_workflow(args)

    if args["which"] == "variant_counting":
        workflow = variant_counting_workflow(args)

    if args["which"] == "ltm":
        workflow = ltm_workflow(args)

    if args["which"] == "multi_sample_pseudo_bulk":
        workflow = multi_sample_workflow(args)

    if workflow:
        pyp.run(workflow)


if __name__ == "__main__":
    main()
