import pypeliner
from cmdline import parse_args
from align import align_workflow
from hmmcopy import hmmcopy_workflow
from aneufinder import aneufinder_workflow
from merge_bams import merge_bams_workflow
from variant_calling import variant_calling_workflow
from germline_calling import germline_calling_workflow
from breakpoint_calling import breakpoint_calling_workflow
from split_bam import split_bam_workflow
from generate_config import generate_config
from clean_sentinels import clean_sentinels
from copy_number import copy_number_calling_workflow

from single_cell.utils import helpers

# from copyclone import copyclone_workflow

def main():

    args = parse_args()

    args = helpers.generate_configs_in_temp(args)

    if "generate_config" in args:
        generate_config(args["generate_config"])

    if "clean_sentinels" in args:
        clean_sentinels(args["clean_sentinels"])

    if "align" in args:
        pyp = pypeliner.app.Pypeline(config=args["align"])
        workflow = pypeliner.workflow.Workflow()
        workflow = align_workflow(workflow, args["align"])
        pyp.run(workflow)

    if "hmmcopy" in args:
        pyp = pypeliner.app.Pypeline(config=args["hmmcopy"])
        workflow = pypeliner.workflow.Workflow()
        workflow = hmmcopy_workflow(workflow, args["hmmcopy"])
        pyp.run(workflow)

#     if "copyclone" in args["modes"]:
#         pyp = pypeliner.app.Pypeline(config=args["copyclone"])
#         workflow = pypeliner.workflow.Workflow()
#         workflow = copyclone_workflow(workflow, args["copyclone"])
#         pyp.run(workflow)

    if "aneufinder" in args:
        pyp = pypeliner.app.Pypeline(config=args["aneufinder"])
        workflow = pypeliner.workflow.Workflow()
        workflow = aneufinder_workflow(workflow, args["aneufinder"])
        pyp.run(workflow)
        
    if "merge_bams" in args:
        pyp = pypeliner.app.Pypeline(config=args["merge_bams"])
        workflow = pypeliner.workflow.Workflow()
        workflow = merge_bams_workflow(workflow, args["merge_bams"])
        pyp.run(workflow)

    if "split_bam" in args:
        pyp = pypeliner.app.Pypeline(config=args["split_bam"])
        workflow = pypeliner.workflow.Workflow()
        workflow = split_bam_workflow(workflow, args["split_bam"])
        pyp.run(workflow)

    if "variant_calling" in args:
        pyp = pypeliner.app.Pypeline(config=args["variant_calling"])
        workflow = pypeliner.workflow.Workflow()
        workflow = variant_calling_workflow(workflow, args["variant_calling"])
        pyp.run(workflow)

    if "copy_number_calling" in args:
        pyp = pypeliner.app.Pypeline(config=args["copy_number_calling"])
        workflow = pypeliner.workflow.Workflow()
        workflow = copy_number_calling_workflow(workflow, args["copy_number_calling"])
        pyp.run(workflow)

    if "germline_calling" in args:
        pyp = pypeliner.app.Pypeline(config=args["germline_calling"])
        workflow = pypeliner.workflow.Workflow()
        workflow = germline_calling_workflow(workflow, args["germline_calling"])
        pyp.run(workflow)

    if "breakpoint_calling" in args:
        pyp = pypeliner.app.Pypeline(config=args["breakpoint_calling"])
        workflow = pypeliner.workflow.Workflow()
        workflow = breakpoint_calling_workflow(workflow, args["breakpoint_calling"])
        pyp.run(workflow)

if __name__ == "__main__":
    main()
