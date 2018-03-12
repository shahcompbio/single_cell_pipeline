import os
import pypeliner
from cmdline import parse_args
from align import align_workflow
from hmmcopy import hmmcopy_workflow
from aneufinder import aneufinder_workflow
from summary import summary_workflow
from pseudo_wgs import pseudo_wgs_workflow
from variant_calling import variant_calling_workflow

def main():

    args = parse_args()

    pyp = pypeliner.app.Pypeline(config=args)

    if args["which"] == "align" or args["which"]=="all" or args["which"]=="qc":
        workflow = pypeliner.workflow.Workflow()
        workflow = align_workflow(workflow, args)
        pyp.run(workflow)

    if args["which"] == "hmmcopy" or args["which"]=="all" or args["which"]=="qc":
        workflow = pypeliner.workflow.Workflow()
        workflow = hmmcopy_workflow(workflow, args)
        pyp.run(workflow)

    if args["which"] == 'aneufinder' or args["which"]=="all":
        workflow = pypeliner.workflow.Workflow()
        workflow = aneufinder_workflow(workflow, args)
        pyp.run(workflow)
    
    if args["which"] == "summary" or args["which"]=="all" or args["which"]=="qc":
        workflow = pypeliner.workflow.Workflow()
        workflow = summary_workflow(workflow, args)
        pyp.run(workflow)
        
    if args["which"] == "pseudo_wgs" or args["which"]=="all":
        workflow = pypeliner.workflow.Workflow()
        workflow = pseudo_wgs_workflow(workflow, args)
        pyp.run(workflow)

    if args["which"] == "variant_calling" or args["which"]=="all":
        workflow = pypeliner.workflow.Workflow()
        workflow = variant_calling_workflow(workflow, args)
        pyp.run(workflow)



if __name__ == '__main__':
    main()
