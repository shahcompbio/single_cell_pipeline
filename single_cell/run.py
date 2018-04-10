import os
import pypeliner
from cmdline import parse_args
from align import align_workflow
from hmmcopy import hmmcopy_workflow
from aneufinder import aneufinder_workflow
from pseudo_wgs import pseudo_wgs_workflow
from variant_calling import variant_calling_workflow
from germline_calling import germline_calling_workflow
from split_normal import split_normal_workflow

# from copyclone import copyclone_workflow

def main():

    args = parse_args()

    pyp = pypeliner.app.Pypeline(config=args)

    if "align" in args["modes"]:
        workflow = pypeliner.workflow.Workflow()
        workflow = align_workflow(workflow, args)
        pyp.run(workflow)

    if "hmmcopy" in args["modes"]:
        workflow = pypeliner.workflow.Workflow()
        workflow = hmmcopy_workflow(workflow, args)
        pyp.run(workflow)

#     if "copyclone" in args["modes"]:
#         workflow = pypeliner.workflow.Workflow()
#         workflow = copyclone_workflow(workflow, args)
#         pyp.run(workflow)

    if "aneufinder" in args["modes"]:
        workflow = pypeliner.workflow.Workflow()
        workflow = aneufinder_workflow(workflow, args)
        pyp.run(workflow)
        
    if "pseudo_wgs" in args["modes"]:
        workflow = pypeliner.workflow.Workflow()
        workflow = pseudo_wgs_workflow(workflow, args, args["merged_wgs_template"], args["input_yaml"], mode="tumour")
        pyp.run(workflow)

    if "split_normal" in args["modes"]:
        workflow = pypeliner.workflow.Workflow()

        if args.get("matched_normal", None):
            workflow = split_normal_workflow(workflow, args)
        else:
            workflow = pseudo_wgs_workflow(workflow, args, args["normal_split_template"], args["normal_yaml"], mode="normal")
        pyp.run(workflow)

    if "variant_calling" in args["modes"]:
        workflow = pypeliner.workflow.Workflow()
        workflow = variant_calling_workflow(workflow, args)
        pyp.run(workflow)

    if "germline_calling" in args["modes"]:
        workflow = pypeliner.workflow.Workflow()
        workflow = germline_calling_workflow(workflow, args)
        pyp.run(workflow)



if __name__ == '__main__':
    main()
