'''
Created on Jul 6, 2017

@author: dgrewal
'''
import pypeliner.managed as mgd

import pypeliner


def create_corrupt_tree_workflow(
        hmmcopy_metrics, hmmcopy_reads,
        tree_newick, consensus_newick, phylo_csv, rank_trees,
        filtered_input, pdf_output, library_id, config
):
    ctx = {'docker_image': config['docker']['single_cell_pipeline']}

    corrupt_tree_params = config['corrupt_tree_params']

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.transform(
        name='generate_corrupt_tree_inputs',
        func="single_cell.workflows.corrupt_tree.tasks.preprocess",
        args=(
            mgd.InputFile(hmmcopy_reads),
            mgd.InputFile(hmmcopy_metrics),
            mgd.TempOutputFile('corrupt_tree_input.csv'),
            mgd.TempSpace("corrupt_tree_gen"),
            library_id,
        ),
        kwargs={'docker_image': config['docker']['corrupt_tree']}
    )

    workflow.transform(
        name='straighten_jitter',
        func="single_cell.workflows.corrupt_tree.tasks.func_straighten",
        args=(
            pypeliner.managed.TempInputFile('corrupt_tree_input.csv'),
            pypeliner.managed.TempOutputFile('jitter_straighten.csv'),
            pypeliner.managed.TempSpace("straighten_jitter_temp"),
        ),
        kwargs={
            'docker_image': config['docker']['corrupt_tree'],
            'neighborhood_size': corrupt_tree_params['neighborhood_size']
        }
    )

    workflow.transform(
        name='data_exploration',
        func="single_cell.workflows.corrupt_tree.tasks.data_exploration",
        args=(
            pypeliner.managed.TempInputFile('jitter_straighten.csv'),
            pypeliner.managed.TempOutputFile('data_exploration.pdf'),
        ),
        kwargs={'docker_image': config['docker']['corrupt_tree']}
    )

    workflow.transform(
        name='filter',
        func="single_cell.workflows.corrupt_tree.tasks.func_filter",
        args=(
            pypeliner.managed.TempInputFile('jitter_straighten.csv'),
            pypeliner.managed.OutputFile(filtered_input),
            pypeliner.managed.TempSpace("filtered"),
        ),
        kwargs={
            'docker_image': config['docker']['corrupt_tree'],
            'lower_fraction': corrupt_tree_params['lower_fraction']
        }
    )

    workflow.transform(
        name='inference',
        func="single_cell.workflows.corrupt_tree.tasks.func_inference",
        args=(
            pypeliner.managed.InputFile(filtered_input),
            pypeliner.managed.OutputFile(phylo_csv),
            pypeliner.managed.TempOutputFile('fnr_out.csv'),
            pypeliner.managed.TempOutputFile('fpr_out.csv'),
            pypeliner.managed.TempOutputFile('log_density_out.csv'),
            pypeliner.managed.TempSpace("inference"),
        ),
        kwargs={
            'docker_image': config['docker']['corrupt_tree'],
            'engine_nchains': corrupt_tree_params['engine_nchains'],
            'engine_nscans': corrupt_tree_params['engine_nscans'],
            'model_fpr_bound': corrupt_tree_params['model_fpr_bound'],
            'model_fnr_bound': corrupt_tree_params['model_fnr_bound'],
        }
    )

    workflow.transform(
        name='find_consensus',
        func="single_cell.workflows.corrupt_tree.tasks.func_find_consensus",
        args=(
            pypeliner.managed.InputFile(phylo_csv),
            pypeliner.managed.OutputFile(consensus_newick),
            pypeliner.managed.TempSpace("consensus"),
        ),
        kwargs={'docker_image': config['docker']['corrupt_tree']}
    )

    workflow.transform(
        name='error_rates_viz_fnr',
        func="single_cell.workflows.corrupt_tree.tasks.error_rates_viz",
        args=(
            pypeliner.managed.TempInputFile('fnr_out.csv'),
            pypeliner.managed.TempOutputFile('trace_fnr.pdf'),
            pypeliner.managed.TempOutputFile('box_fnr.pdf'),
        ),
        kwargs={'docker_image': config['docker']['corrupt_tree']}
    )

    workflow.transform(
        name='error_rates_viz_fpr',
        func="single_cell.workflows.corrupt_tree.tasks.error_rates_viz",
        args=(
            pypeliner.managed.TempInputFile('fpr_out.csv'),
            pypeliner.managed.TempOutputFile('trace_fpr.pdf'),
            pypeliner.managed.TempOutputFile('box_fpr.pdf'),
        ),
        kwargs={'docker_image': config['docker']['corrupt_tree']}
    )

    workflow.transform(
        name='error_rates_viz_logd',
        func="single_cell.workflows.corrupt_tree.tasks.logd",
        args=(
            pypeliner.managed.TempInputFile('log_density_out.csv'),
            pypeliner.managed.TempOutputFile('trace_logd.pdf'),
        ),
        kwargs={'docker_image': config['docker']['corrupt_tree']}
    )

    workflow.transform(
        name='average_tip_indicators',
        func="single_cell.workflows.corrupt_tree.tasks.func_average",
        args=(
            pypeliner.managed.InputFile(phylo_csv),
            pypeliner.managed.TempOutputFile('average_tip.csv'),
            pypeliner.managed.TempSpace("average"),
        ),
        kwargs={'docker_image': config['docker']['corrupt_tree']}
    )

    workflow.transform(
        name='viz_edge_locus_decode',
        func="single_cell.workflows.corrupt_tree.tasks.func_viz_edge_locus_decode",
        args=(
            pypeliner.managed.TempInputFile('average_tip.csv'),
            pypeliner.managed.InputFile(filtered_input),
            pypeliner.managed.InputFile(consensus_newick),
            pypeliner.managed.TempOutputFile('viz_edge_locus.pdf'),
            pypeliner.managed.TempSpace("viz_edge_locus_decode"),
        ),
        kwargs={'docker_image': config['docker']['corrupt_tree']}
    )

    workflow.transform(
        name='decode',
        func="single_cell.workflows.corrupt_tree.tasks.func_decode",
        args=(
            pypeliner.managed.TempInputFile('average_tip.csv'),
            pypeliner.managed.OutputFile(tree_newick),
            pypeliner.managed.TempOutputFile('tree_decode_out.csv'),
            pypeliner.managed.TempSpace("decode"),
        ),
        kwargs={'docker_image': config['docker']['corrupt_tree']}
    )

    workflow.transform(
        name='viz',
        func="single_cell.workflows.corrupt_tree.tasks.func_viz",
        args=(
            pypeliner.managed.TempInputFile('average_tip.csv'),
            pypeliner.managed.InputFile(filtered_input),
            pypeliner.managed.InputFile(tree_newick),
            pypeliner.managed.TempOutputFile("func_viz_out.pdf"),
            pypeliner.managed.TempSpace("viz"),
        ),
        kwargs={'docker_image': config['docker']['corrupt_tree']}
    )

    workflow.transform(
        name='rank_loci',
        func="single_cell.workflows.corrupt_tree.tasks.func_rank_loci",
        args=(
            pypeliner.managed.InputFile(tree_newick),
            pypeliner.managed.InputFile(filtered_input),
            pypeliner.managed.OutputFile(rank_trees),
            pypeliner.managed.TempSpace("rank_loci"),
        ),
        kwargs={'docker_image': config['docker']['corrupt_tree']}
    )

    workflow.transform(
        name='detailed_viz',
        func="single_cell.workflows.corrupt_tree.tasks.func_detailed_viz",
        args=(
            pypeliner.managed.InputFile(rank_trees),
            pypeliner.managed.TempInputFile('average_tip.csv'),
            pypeliner.managed.InputFile(filtered_input),
            pypeliner.managed.TempOutputFile("detailed.pdf"),
            pypeliner.managed.TempSpace("temp_detailed_viz"),
        ),
        kwargs={'docker_image': config['docker']['corrupt_tree']}
    )

    workflow.transform(
        name='merge_all_pdfs',
        func='single_cell.utils.pdfutils.merge_pdfs',
        args=(
            [
                pypeliner.managed.TempInputFile("func_viz_out.pdf"),
                pypeliner.managed.TempInputFile('trace_fnr.pdf'),
                pypeliner.managed.TempInputFile('box_fnr.pdf'),
                pypeliner.managed.TempInputFile('trace_fpr.pdf'),
                pypeliner.managed.TempInputFile('box_fpr.pdf'),
                pypeliner.managed.TempInputFile("detailed.pdf"),
                pypeliner.managed.TempInputFile('trace_logd.pdf'),
                pypeliner.managed.TempInputFile('data_exploration.pdf'),
                pypeliner.managed.TempInputFile('viz_edge_locus.pdf'),
            ],
            pypeliner.managed.OutputFile(pdf_output)
        )
    )

    return workflow
