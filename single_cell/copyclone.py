'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import pypeliner.managed as mgd
from workflows import copyclone
from utils import helpers


def parameter_set1():

    A = [0.994, 0.994, 0.994, 0.994, 0.994, 0.994, 0.994]
    alpha_A = [1000, 1000, 1000, 1000, 1000, 1000, 1000]
    pi = [0.05, 0.1, 0.5, 0.2, 0.05, 0.05, 0.05]
    alpha_pi = [2, 2, 50, 2, 2, 2, 2]
    tau = [500, 25, 25, 25, 25, 25, 15]
    nu = [5, 5, 5, 5, 5, 5, 5]
    eta = [5000, 5000, 5000, 5000, 5000, 5000, 5000]
    shape = [3, 30, 30, 30, 30, 30, 20]
    rate = [0.01, 1, 1, 1, 1, 1, 1]
    ploidy_states = [2, 3, 4]

    num_states = 7

    params = {'A': A, 'alpha_A': alpha_A, 'pi': pi, 'alpha_pi': alpha_pi,
              'tau': tau, 'nu': nu, 'eta': eta, 'shape': shape, 'rate': rate,
              'ploidy_states': ploidy_states, 'num_states': num_states}

    return ('set1', params)


def parameter_set2():

    A = [0.994, 0.994, 0.994, 0.9994, 0.994, 0.994, 0.994]
    alpha_A = [1000, 1000, 1000, 10000, 1000, 1000, 1000]
    pi = [0.05, 0.1, 0.5, 0.2, 0.05, 0.05, 0.05]
    alpha_pi = [2, 2, 50, 2, 2, 2, 2]
    tau = [500, 25, 25, 25, 25, 25, 15]
    nu = [5, 5, 5, 5, 5, 5, 5]
    eta = [5000, 5000, 5000, 5000, 5000, 5000, 5000]
    shape = [3, 30, 30, 30, 30, 30, 20]
    rate = [0.01, 1, 1, 1, 1, 1, 1]
    ploidy_states = [2, 3, 4]

    num_states = 7

    params = {'A': A, 'alpha_A': alpha_A, 'pi': pi, 'alpha_pi': alpha_pi,
              'tau': tau, 'nu': nu, 'eta': eta, 'shape': shape, 'rate': rate,
              'ploidy_states': ploidy_states, 'num_states': num_states}

    return ('set2', params)


def copyclone_workflow(workflow, args):

    config = helpers.load_config(args)
    sampleids = helpers.get_samples(args['input_yaml'])
    bam_files, bai_files = helpers.get_bams(args['input_yaml'])

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sampleids,
    )

    libid = args['library_id']

    parameters = [parameter_set1(), parameter_set2()]

    h5_files = []

    optimal = os.path.join(
        args['out_dir'],
        'results',
        "copyclone",
        "optimal.txt")

    for name, param in parameters:

        results_dir = os.path.join(
            args['out_dir'],
            'results',
            "copyclone",
            name)
        output_filename = os.path.join(
            results_dir,
            '{}_copyclone.h5'.format(libid))

        h5_files.append(output_filename)

        plots_dir = os.path.join(results_dir, "plots")
        segs_pdf = os.path.join(plots_dir, '{}_segments.pdf'.format(libid))
        reads_pdf = os.path.join(plots_dir, '{}_bias.pdf'.format(libid))
        metrics_pdf = os.path.join(plots_dir, '{}_metrics.pdf'.format(libid))
        kde_pdf = os.path.join(
            plots_dir,
            '{}_kernel_density.pdf'.format(libid))
        heatmap_pdf = os.path.join(plots_dir, '{}_heatmap.pdf'.format(libid))
        heatmap_filt_pdf = os.path.join(
            plots_dir,
            '{}_heatmap_filtered.pdf'.format(libid))

        # monkeypatch params
        config['copyclone']['A'] = param['A']
        config['copyclone']['alpha_A'] = param['alpha_A']
        config['copyclone']['pi'] = param['pi']
        config['copyclone']['alpha_pi'] = param['alpha_pi']
        config['copyclone']['tau'] = param['tau']
        config['copyclone']['nu'] = param['nu']
        config['copyclone']['eta'] = param['eta']
        config['copyclone']['shape'] = param['shape']
        config['copyclone']['rate'] = param['rate']
        config['copyclone']['ploidy_states'] = param['ploidy_states']
        config['copyclone']['num_states'] = param['num_states']

        workflow.subworkflow(
            name='copyclone_workflow_' + name,
            func=copyclone.create_copyclone_workflow,
            args=(mgd.InputFile('bam_markdups', 'sample_id', fnames=bam_files),
                  mgd.InputFile(
                'bam_markdups_index',
                'sample_id',
                fnames=bai_files),
                mgd.OutputFile(output_filename),
                mgd.OutputFile(segs_pdf),
                mgd.OutputFile(reads_pdf),
                mgd.OutputFile(metrics_pdf),
                mgd.OutputFile(kde_pdf),
                mgd.OutputFile(heatmap_pdf),
                mgd.OutputFile(heatmap_filt_pdf),
                sampleids,
                config,
                args,
                results_dir
            ),
        )

    workflow.transform(
        name='getoptimal',
        ctx={
            'mem': config["memory"]['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1},
        func=compute_optimal,
        args=(
            [mgd.InputFile(h5file) for h5file in h5_files],
            parameters,
            mgd.OutputFile(optimal)
        ),
    )

    return workflow


def compute_optimal(h5files, parameters, optimal):
    import pandas as pd
    import numpy as np

    lowest_halfiness = None
    for h5data, (name, _) in zip(h5files, parameters):
        metrics = pd.read_hdf(h5data, '/copyclone/metrics')

        halfiness = np.nanmean(metrics.scaled_halfiness)

        if not lowest_halfiness:
            lowest_halfiness = (name, halfiness)

        if lowest_halfiness[1] > halfiness:
            lowest_halfiness = (name, halfiness)

        print halfiness

    with open(optimal, "w") as output:
        output.write("\t".join(map(str, lowest_halfiness)))
