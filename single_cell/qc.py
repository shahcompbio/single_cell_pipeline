import os
import re

import pypeliner.managed as mgd
from single_cell.utils import inpututils
from single_cell.workflows import align

import pypeliner
import sys




def qc_workflow(args):
    yaml = inpututils.load_yaml(args["input_yaml"])
    config = inpututils.load_config(args)
    
    config = config["sv_genotyping"]


    #inputs
    jira_ids =  {k: v["jira_id"] for k, v in yaml.items()}
    snv_jira_ids = {k: v["snv_genotyping_jira_id"]  for k,v in yaml.items()}
    library_ids = {k: v["library_id"] for k, v in yaml.items()}
    cell_ids = {k: k for k, v in yaml.items()}
    # pseudobulk_groups = {k: cell + "_" + lib for (k, cell), (k,lib) in zip(cell_ids, library_ids) }
    out_dir = args["out_dir"]
    tmp_dir = args["tmpdir"]
    #outputs
    reporthtmls = {k: os.path.join(out_dir, k + "_main_report.html") for k, v in yaml.items()}
    mutationreporthtmls = {k: os.path.join(out_dir, k + "_mutation_report.html") for k, v in yaml.items()}
    mafs = {k: os.path.join(out_dir, k + "_maf.maf") for k, v in yaml.items()}



    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config['docker']['single_cell_pipeline']}
    )

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=list(yaml.keys()),
    )

    workflow.subworkflow(
        name='create_qc_workflow',
        ctx={'mem':450, 'docker_image': config['docker']['single_cell_pipeline']},
        func="single_cell.workflows.qc.create_qc_workflow",
        axes=('cell_id',),
        args=(
            mgd.InputFile('cell_id', 'cell_id', fnames=cell_ids),
            mgd.InputFile('jira_id', 'cell_id', fnames=jira_ids),
            mgd.InputFile('snv_jira_id', 'cell_id', fnames=snv_jira_ids),
            mgd.InputFile('library_id', 'cell_id', fnames=library_ids),
            tmp_dir,
            out_dir,
            mgd.InputFile('reporthtmls', 'cell_id', fnames=reporthtmls),
            mgd.InputFile('mutationreporthtmls', 'cell_id', fnames=mutationreporthtmls),
            mgd.InputFile('mafs', 'cell_id', fnames=mafs),
        ),
    )
    # workflow.subworkflow(
    #     name='create_qc_workflow',
    #     ctx={'mem':450, 'docker_image': config['docker']['single_cell_pipeline']},
    #     func="single_cell.workflows.qc.create_qc_workflow",
    #     args=(
    #         cell_ids,
    #         library_ids,
    #         snv_jira_ids,
    #         library_ids,
    #         pseudobulk_groups,xw
    #         pseudobulk_groups,
    #         tmp_dir,
    #         out_dir,
    #         reporthtmls,
    #         mutationreporthtmls,
    #         mafs,
    #     ),
    # )

    

    # workflow.transform(
    #     name='generate_meta_files_results',
    #     func='wgs.utils.helpers.generate_and_upload_metadata',
    #     args=(
    #         sys.argv[0:],
    #         args["out_dir"],
    #         outputted_filenames,
    #         mgd.OutputFile(meta_yaml)
    #     ),
    #     kwargs={
    #         'input_yaml_data': helpers.load_yaml(args['input_yaml']),
    #         'input_yaml': mgd.OutputFile(input_yaml_blob),
    #         'metadata': {'type': 'postprocessing'}
    #     }
    # )


    # normals = {sample: yamldata[sample]['normal_bam'] for sample in samples}
    # tumours = {sample: yamldata[sample]['tumour_bam'] for sample in samples}

    # titan = {sample: yamldata[sample]['titan'] for sample in samples}
    # remixt = {sample: yamldata[sample]['remixt'] for sample in samples}
    # breakpoints_consensus = {sample: yamldata[sample]['breakpoints_consensus'] for sample in samples}
    # roh = {sample: yamldata[sample]['roh'] for sample in samples}
    # germline_calls = {sample: yamldata[sample]['germline_calls'] for sample in samples}
    # somatic_calls = {sample: yamldata[sample]['somatic_calls'] for sample in samples}

    # out_dir = args['out_dir']

    # meta_yaml = os.path.join(out_dir, 'pipeline_metadata.yaml')
    # input_yaml_blob = os.path.join(out_dir, 'input.yaml')

    # circos_plot_remixt = os.path.join(out_dir, '{sample_id}', '{sample_id}_circos_remixt.pdf')
    # circos_plot_titan = os.path.join(out_dir, '{sample_id}', '{sample_id}_circos_titan.pdf')

    # genome_wide_plot = os.path.join(out_dir, '{sample_id}', '{sample_id}_genome_wide.pdf')

    # pyp = pypeliner.app.Pypeline(config=args)
    # workflow = pypeliner.workflow.Workflow(
    #     ctx=helpers.get_default_ctx(docker_image=config.containers('wgs'))
    # )

    # workflow.setobj(
    #     obj=mgd.OutputChunks('sample_id'),
    #     value=samples
    # )

    # workflow.subworkflow(
    #     name="postprocessing",
    #     func=postprocessing.create_postprocessing_workflow,
    #     ctx=helpers.get_default_ctx(),
    #     axes=('sample_id',),
    #     args=(
    #         mgd.InputFile('normal.bam', 'sample_id', fnames=normals),
    #         mgd.InputFile('tumour.bam', 'sample_id', fnames=tumours),
    #         titan,
    #         remixt,
    #         breakpoints_consensus,
    #         roh,
    #         germline_calls,
    #         somatic_calls,
    #         mgd.OutputFile('circos_plot_remixt.pdf', 'sample_id', template=circos_plot_remixt),
    #         mgd.OutputFile('circos_plot_titan.pdf', 'sample_id', template=circos_plot_titan),
    #         mgd.OutputFile('genome_wide_plot.pdf', 'sample_id', template=genome_wide_plot),
    #         args['refdir'],
    #         mgd.InputInstance('sample_id'),
    #     ),
    #     kwargs={'single_node': args['single_node']}
    # )

    # outputted_filenames = helpers.expand_list([circos_plot_remixt, circos_plot_titan, genome_wide_plot], samples, "sample_id")

    # workflow.transform(
    #     name='generate_meta_files_results',
    #     func='wgs.utils.helpers.generate_and_upload_metadata',
    #     args=(
    #         sys.argv[0:],
    #         args["out_dir"],
    #         outputted_filenames,
    #         mgd.OutputFile(meta_yaml)
    #     ),
    #     kwargs={
    #         'input_yaml_data': helpers.load_yaml(args['input_yaml']),
    #         'input_yaml': mgd.OutputFile(input_yaml_blob),
    #         'metadata': {'type': 'postprocessing'}
    #     }
    # )

    # pyp.run(workflow)

    # def qc_pipeline(args):
    #     pyp = pypeliner.app.Pypeline(config=args)

    #     workflow = qc_workflow(args)

    #     pyp.run(workflow)
    return workflow
def qc_pipeline(args):
    pyp = pypeliner.app.Pypeline(config=args)

    workflow = qc_workflow(args)

    pyp.run(workflow)