import os
import sys
import argparse

import pypeliner
import pypeliner.managed as mgd


if __name__ == '__main__':

    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    pypeliner.app.add_arguments(argparser)

    parser.add_argument('config_file',
                        help='''Path to yaml config file.''')

    args = vars(argparser.parse_args())

    config = {}
    if args['config'] is not None:
        execfile(args['config'], {}, config)

    pyp = pypeliner.app.Pypeline(config=args)

    #=======================================================================================================================
    # Read Input Files
    #=======================================================================================================================
    try:
        with open(args.config_file) as file:
            config = yaml.load(file)
            
    except IOError as e:
        print 'Unable to open config file: {0}'.format(args.config_file)

    sample_sheet_file = os.path.join(config['nextseq_dir'], 'SampleSheet.csv')

    try:
        with open(sample_sheet_file) as file:
            lines = [x.strip('\n').strip(',') for x in file.readlines()]
        
        run_id = [s.split(',')[1] for s in lines if 'Experiment Name,' in s][0]
        
        library_id = [s.split(',')[1] for s in lines if 'Description,' in s][0]
        
    except IOError as e:
        print 'Unable to open file \'SampleSheet.csv\' in directory: {0}'.format(config['nextseq_dir'])

    workflow = pypeliner.workflow.Workflow()

    







    experiment_template = os.path.join(args['raw_data_dir'], '{sim_id}', 'experiment.pickle')
    experiment_plot_template = os.path.join(args['raw_data_dir'], '{sim_id}', 'experiment_plot.pdf')
    results_template = os.path.join(args['raw_data_dir'], '{sim_id}', 'results.h5')
    evaluation_template = os.path.join(args['raw_data_dir'], '{sim_id}', 'evaluation.h5')

    workflow = pypeliner.workflow.Workflow(default_ctx={'mem': 4})

    workflow.transform(
        name='read_sim_defs',
        ctx={'local': True},
        func=remixt.simulations.pipeline.create_simulations,
        ret=mgd.TempOutputObj('sim_defs', 'sim_id'),
        args=(
            mgd.InputFile(args['sim_defs']),
            config,
            args['ref_data_dir'],
        ),
    )

    workflow.transform(
        name='simulate_experiment',
        axes=('sim_id',),
        func=remixt.simulations.pipeline.simulate_experiment,
        args=(
            mgd.OutputFile('experiment', 'sim_id', template=experiment_template),
            mgd.OutputFile('experiment_plot', 'sim_id', template=experiment_plot_template),
            mgd.TempInputObj('sim_defs', 'sim_id'),
        ),
    )

    if args['simulate_only']:
        pyp.run(workflow)
        sys.exit()

    workflow.subworkflow(
        name='create_tool_workflow',
        axes=('sim_id',),
        func=remixt.workflow.create_fit_model_workflow,
        args=(
            mgd.InputFile('experiment', 'sim_id', template=experiment_template),
            mgd.OutputFile('results', 'sim_id', template=results_template),
            config,
            args['ref_data_dir'],
        ),
    )

    workflow.transform(
        name='evaluate_results',
        axes=('sim_id',),
        func=remixt.simulations.pipeline.evaluate_results_task,
        args=(
            mgd.OutputFile('evaluation', 'sim_id', template=evaluation_template),
            mgd.InputFile('results', 'sim_id', template=results_template),
        ),
        kwargs={
            'experiment_filename': mgd.InputFile('experiment', 'sim_id', template=experiment_template),
        },
    )

    workflow.transform(
        name='merge_evaluations',
        func=remixt.simulations.pipeline.merge_evaluations,
        args=(
            mgd.OutputFile(args['table']),
            mgd.TempInputObj('sim_defs', 'sim_id'),
            mgd.InputFile('evaluation', 'sim_id', template=evaluation_template),
            ['sim_id'],
        ),
    )

    pyp.run(workflow)

