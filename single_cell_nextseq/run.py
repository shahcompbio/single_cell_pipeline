import os
import argparse
import utils
import pypeliner
import pypeliner.managed as mgd
from single_cell_nextseq.workflows import tasks
from workflows import alignment_workflow, hmmcopy_workflow
from workflows import summary_workflow, fastqc_workflow



scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
plot_hmmcopy_script = os.path.join(scripts_directory, 'plot_hmmcopy.py')



def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    pypeliner.app.add_arguments(parser)

    parser.add_argument('nextseq_dir',
                        help='''Path to input nextseq directory.''')

    parser.add_argument('samplesheet',
                        help='''Path to input nextseq directory.''')

    parser.add_argument('out_dir',
                        help='''Path to output files.''')

    parser.add_argument('config_file',
                        help='''Path to yaml config file.''')

    args = vars(parser.parse_args())
    args['tmpdir'] = os.path.join(args['out_dir'], 'tmp')

    return args


def main():

    args = parse_args()

    pyp = pypeliner.app.Pypeline(config=args)

    config = utils.load_config(args)

    fastq_directory = os.path.join(args['out_dir'], 'fastq')


    run_id, library_id, sample_ids, fastq_1_filenames, fastq_2_filenames = utils.read_samplesheet(args,fastq_directory)

    rg_template = utils.get_readgroup_template(library_id, run_id, config)

    trimgalore_results_template_r1 = os.path.join(args['out_dir'], 'trim', '{sample_id}_R1.fastq.gz')
    trimgalore_results_template_r2 = os.path.join(args['out_dir'], 'trim', '{sample_id}_R2.fastq.gz')
    metrics_directory = os.path.join(args['out_dir'], 'metrics')
    metrics_summary_template = os.path.join(metrics_directory, 'summary', '{sample_id}_summary.csv')
    metrics_gcmatrix_template = os.path.join(metrics_directory, 'summary', '{sample_id}_gcmatrix.csv')
    cnmatrix_template = os.path.join(metrics_directory, 'summary', '{sample_id}_cnmatrix.csv')


    bam_directory = os.path.join(args['out_dir'], 'bams')
    bam_template = os.path.join(bam_directory, '{sample_id}.bam')
    bam_index_template = os.path.join(bam_directory, '{sample_id}.bam.bai')

    hmmcopy_directory = os.path.join(args['out_dir'], 'hmmcopy')
    hmmcopy_wig_template = os.path.join(hmmcopy_directory, '{sample_id}_readcount.wig')
    hmmcopy_reads_template = os.path.join(hmmcopy_directory, '{sample_id}_reads.csv')
    hmmcopy_segments_template = os.path.join(hmmcopy_directory, '{sample_id}_segments.csv')
    hmmcopy_parameters_template = os.path.join(hmmcopy_directory, '{sample_id}_parameters.csv')
#     hmmcopy_posteriors_template = os.path.join(hmmcopy_directory, '{sample_id}_posteriors.csv')
    hmmcopy_hmm_metrics_template = os.path.join(hmmcopy_directory, '{sample_id}_hmm_metrics.csv')


    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )

    workflow.transform(
        name='demultiplex_fastq_files',
        ctx={'mem': 64},
        func=single_cell_nextseq.tasks.demultiplex_fastq_files,
        args=(
            mgd.InputFile(args['samplesheet']),
            args['nextseq_dir'],
            mgd.OutputFile('fastq_1', 'sample_id', fnames=fastq_1_filenames, axes_origin=[]),
            mgd.OutputFile('fastq_2', 'sample_id', fnames=fastq_2_filenames, axes_origin=[]),
            mgd.TempSpace('demultiplex_temp'),
            config['bcl2fastq']
        ),
    )


    workflow.subworkflow(
        name='fastqc_workflow',
        axes=('sample_id',),
        func=fastqc_workflow.create_fastqc_workflow,
        args=(mgd.InputFile('fastq_1', 'sample_id', fnames=fastq_1_filenames),
              mgd.InputFile('fastq_2', 'sample_id', fnames=fastq_2_filenames),
              mgd.OutputFile('fastq_trim_1', 'sample_id', template=trimgalore_results_template_r1),
              mgd.OutputFile('fastq_trim_2', 'sample_id', template=trimgalore_results_template_r2),
              config,
              metrics_directory,
              mgd.InputInstance('sample_id'),
            ),
        )


    workflow.subworkflow(
        name='alignment_workflow',
        axes=('sample_id',),
        func=alignment_workflow.create_alignment_workflow,
        args=(
            mgd.InputFile('fastq_trim_1', 'sample_id', template=trimgalore_results_template_r1),
            mgd.InputFile('fastq_trim_2', 'sample_id', template=trimgalore_results_template_r2),
            mgd.OutputFile('bam', 'sample_id', template=bam_template),
            mgd.OutputFile('bam_index', 'sample_id', template=bam_index_template),
            mgd.InputFile(config['ref_genome']),
            mgd.Template(rg_template, 'sample_id'),
            mgd.OutputFile(metrics_summary_template, 'sample_id'),
            mgd.OutputFile(metrics_gcmatrix_template, 'sample_id'),            
            metrics_directory,
            mgd.InputInstance('sample_id'),
            mgd.InputFile(args['samplesheet']),
            config,
        ),
    )

    workflow.subworkflow(
        name='hmmcopy_workflow',
        axes=('sample_id',),
        func=hmmcopy_workflow.create_hmmcopy_workflow,
        args=(
            mgd.InputFile('bam', 'sample_id', template=bam_template),
            mgd.OutputFile(hmmcopy_wig_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_reads_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_segments_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_parameters_template, 'sample_id'),
#             mgd.OutputFile(hmmcopy_posteriors_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_hmm_metrics_template, 'sample_id'),
            mgd.OutputFile(cnmatrix_template, 'sample_id'),
            mgd.InputInstance('sample_id'),
            config,
            args
        ),
    )


    workflow.subworkflow(
        name='summary_workflow',
        func=summary_workflow.create_summary_workflow,
        args=(
            mgd.InputFile(hmmcopy_segments_template, 'sample_id'),
            mgd.InputFile(hmmcopy_reads_template, 'sample_id'),
            mgd.InputFile(hmmcopy_hmm_metrics_template, 'sample_id'),
            mgd.InputFile(metrics_summary_template, 'sample_id'),
            mgd.InputFile(metrics_gcmatrix_template, 'sample_id'),
            mgd.InputFile(cnmatrix_template, 'sample_id'),
            config,
            args,
        ),
    )

    pyp.run(workflow)

if __name__ == '__main__':
    main()
