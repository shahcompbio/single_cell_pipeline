'''
Created on Jul 24, 2017

@author: dgrewal
'''
import glob
import os
import shutil

from pypeliner.commandline import execute
from single_cell.utils import helpers
from single_cell.utils import pdfutils


def run_in_dir(command, outputs, cmd_outputs, dirpath, docker_image=None):
    outputs = [os.path.abspath(val) for val in outputs]

    og_dir = os.getcwd()

    helpers.makedirs(dirpath)
    os.chdir(dirpath)

    execute(*command, docker_image=docker_image)

    for local_path, final_path in zip(cmd_outputs, outputs):
        shutil.copy(local_path, final_path)

    os.chdir(og_dir)


def preprocess(reads, metrics, output, tempdir, tag, docker_image=None):
    execute(
        'corrupt_tree_preprocessing.R',
        '--reads', reads,
        '--metrics', metrics,
        '--outdir', tempdir,
        '--datatag', tag,
        docker_image=docker_image
    )

    outfile = os.path.join(tempdir, tag, '{}_cnvs_corrupt_no_padding.csv'.format(tag.lower()))

    shutil.copyfile(outfile, output)


def data_exploration(straighten_jitter, output, docker_image=None):
    execute(
        'data_exploration.R',
        '-f',
        straighten_jitter,
        '-o',
        output,
        docker_image=docker_image
    )


def func_straighten(infile, outfile, tmpdir, docker_image=None, neighborhood_size=2):
    cmd = ('corrupt-straighten',
           '--experimentConfigs.saveStandardStreams',
           'false',
           '--experimentConfigs.recordGitInfo',
           'false',
           '--experimentConfigs.managedExecutionFolder',
           'false',
           '--input',
           os.path.abspath(infile),
           '--neighborhoodSize',
           neighborhood_size,
           )

    run_in_dir(cmd, [outfile], ['output.csv'], tmpdir, docker_image=docker_image)


def func_filter(infile, outfile, tmpdir, docker_image=None, lower_fraction=0.05):
    infile = os.path.abspath(infile)

    cmd = ('corrupt-filter',
           '--experimentConfigs.saveStandardStreams',
           'false',
           '--experimentConfigs.recordGitInfo',
           'false',
           '--experimentConfigs.managedExecutionFolder',
           'false',
           '--input',
           infile,
           '--lowerFraction',
           lower_fraction,
           )

    run_in_dir(cmd, [outfile], ['filtered.csv'], tmpdir, docker_image=docker_image)


def func_inference(
        infile, phylo_outfile, fnr_outfile, fpr_outfile,
        log_density_outfile, tmpdir,
        docker_image=None,
        engine_nchains=1,
        engine_nscans=10000,
        model_fpr_bound=0.1,
        model_fnr_bound=0.5
):
    infile = os.path.abspath(infile)

    cmd = ('corrupt-infer-with-noisy-params',
           '--experimentConfigs.saveStandardStreams',
           'false',
           '--experimentConfigs.recordGitInfo',
           'false',
           '--experimentConfigs.managedExecutionFolder',
           'false',
           '--model.globalParameterization',
           'true',
           '--model.binaryMatrix',
           infile,
           '--model.fprBound',
           model_fpr_bound,
           '--model.fnrBound',
           model_fnr_bound,
           '--engine',
           'PT',
           '--engine.ladder',
           'Polynomial',
           '--engine.nScans',
           engine_nscans,
           '--engine.nPassesPerScan',
           '1',
           '--engine.nChains',
           engine_nchains,
           '--engine.nThreads',
           'Fixed',
           '--engine.nThreads.number',
           '1',
           )

    run_in_dir(
        cmd,
        [phylo_outfile, fnr_outfile, fpr_outfile, log_density_outfile],
        ["./samples/phylo.csv", "./samples/fnr.csv", "./samples/fpr.csv", "./samples/logDensity.csv"],
        tmpdir,
        docker_image=docker_image
    )


def func_find_consensus(infile, outfile, tmpdir, docker_image=None):
    infile = os.path.abspath(infile)

    cmd = ('corrupt-l1-decode',
           '--experimentConfigs.saveStandardStreams',
           'false',
           '--experimentConfigs.recordGitInfo',
           'false',
           '--experimentConfigs.managedExecutionFolder',
           'false',
           '--samples',
           infile,
           )

    run_in_dir(cmd, [outfile], ['consensus.newick'], tmpdir, docker_image=docker_image)


def error_rates_viz(fnr_csv, trace_fnr_pdf, box_fnr_pdf, docker_image=None):
    execute(
        'error_rates_viz.R',
        '-f',
        fnr_csv,
        '-t',
        trace_fnr_pdf,
        '-b',
        box_fnr_pdf,
        docker_image=docker_image
    )


def logd(log_density_csv, trace_log_density_pdf, docker_image=None):
    execute(
        'logd.R',
        '-f',
        log_density_csv,
        '-o',
        trace_log_density_pdf,
        docker_image=docker_image
    )


def func_average(infile, outfile, tmpdir, docker_image=None):
    infile = os.path.abspath(infile)
    cmd = ('corrupt-average',
           '--experimentConfigs.saveStandardStreams',
           'false',
           '--experimentConfigs.recordGitInfo',
           'false',
           '--experimentConfigs.managedExecutionFolder',
           'false',
           '--csvFile',
           infile,
           '--logisticTransform',
           'false',
           )

    run_in_dir(cmd, [outfile], ['average.csv'], tmpdir, docker_image=docker_image)


def func_viz_edge_locus_decode(infile_average, infile_filtered, infile_consensus, outfile, tmpdir, docker_image=None):
    infile_average = os.path.abspath(infile_average)
    infile_consensus = os.path.abspath(infile_consensus)
    infile_filtered = os.path.abspath(infile_filtered)

    cmd = (
        'xvfb-run',
        'corrupt-viz',
        '--experimentConfigs.saveStandardStreams',
        'false',
        '--experimentConfigs.recordGitInfo',
        'false',
        '--experimentConfigs.managedExecutionFolder',
        'false',
        '--matrices',
        infile_average,
        infile_filtered,
        '--phylo',
        'file',
        infile_consensus,
        '--size',
        'width',
        '300',
    )

    run_in_dir(cmd, [outfile], ['output.pdf'], tmpdir, docker_image=docker_image)


def func_decode(infile, outfile_newick, outfile_csv, tmpdir, docker_image=None):
    infile = os.path.abspath(infile)

    cmd = ('corrupt-greedy',
           '--experimentConfigs.saveStandardStreams',
           'false',
           '--experimentConfigs.recordGitInfo',
           'false',
           '--experimentConfigs.managedExecutionFolder',
           'false',
           '--tipInclusionProbabilities',
           'ReadOnlyCLMatrix',
           infile,
           )

    run_in_dir(
        cmd,
        [outfile_newick, outfile_csv],
        ['tree.newick', 'trees.csv'],
        tmpdir,
        docker_image=docker_image)


def func_viz(infile_average, infile_filtered, infile_tree, outfile, tmpdir, docker_image=None):
    infile_filtered = os.path.abspath(infile_filtered)
    infile_average = os.path.abspath(infile_average)
    infile_tree = os.path.abspath(infile_tree)

    cmd = (
        'xvfb-run',
        'corrupt-viz',
        '--experimentConfigs.saveStandardStreams',
        'false',
        '--experimentConfigs.recordGitInfo',
        'false',
        '--experimentConfigs.managedExecutionFolder',
        'false',
        '--matrices',
        infile_average,
        infile_filtered,
        '--phylo',
        'file',
        infile_tree,
        '--size',
        'width',
        '300',
    )

    run_in_dir(cmd, [outfile], ['output.pdf'], tmpdir, docker_image=docker_image)


def func_rank_loci(infile_tree, infile_filtered, outfile, tmpdir, docker_image=None):
    infile_tree = os.path.abspath(infile_tree)
    infile_filtered = os.path.abspath(infile_filtered)

    cmd = ('corrupt-rank-loci',
           '--experimentConfigs.saveStandardStreams',
           'false',
           '--experimentConfigs.recordGitInfo',
           'false',
           '--experimentConfigs.managedExecutionFolder',
           'false',
           '--phylo',
           'file',
           infile_tree,
           '--binaryMatrix',
           infile_filtered,
           '--thinningPeriod',
           '10',
           )

    run_in_dir(cmd, [outfile], ['trees.csv'], tmpdir, docker_image=docker_image)


def func_detailed_viz(infile_rank_tree, infile_average, infile_filtered, detailed_pdf, tempdir, docker_image=None):
    infile_filtered = os.path.abspath(infile_filtered)
    infile_average = os.path.abspath(infile_average)
    infile_rank_tree = os.path.abspath(infile_rank_tree)

    cmd = (
        'xvfb-run',
        'corrupt-viz-growth',
        '--experimentConfigs.saveStandardStreams',
        'false',
        '--experimentConfigs.recordGitInfo',
        'false',
        '--experimentConfigs.managedExecutionFolder',
        'false',
        '--phylogenies',
        infile_rank_tree,
        '--matrices',
        infile_average,
        infile_filtered,
        '--size',
        'width',
        '300',
    )

    run_in_dir(cmd, [], [], tempdir, docker_image=docker_image)

    pdffiles = glob.glob('{}/*.pdf'.format(tempdir))

    pdfutils.merge_pdfs(pdffiles, detailed_pdf)
