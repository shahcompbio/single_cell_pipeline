import os
import sys
import argparse

import pypeliner
import pypeliner.managed as mgd

import single_cell_nextseq.tasks
import single_cell_nextseq.workflow

import utils as utl


scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
plot_hmmcopy_script = os.path.join(scripts_directory, 'plot_hmmcopy.py')


if __name__ == '__main__':
    main()


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

    config = utl.load_config(args)

    fastq_directory = os.path.join(args['out_dir'], 'fastq')

    run_id, library_id, sample_ids, fastq_1_filenames, fastq_2_filenames, fastq_1_basenames, fastq_2_basenames = utl.read_samplesheet(args,fastq_directory)

    rg_template = utl.get_readgroup_template(library_id, run_id, config)

    outdir = args['out_dir']


    sample_info_filename = os.path.join(args['out_dir'], 'sample_info.csv')
    trimgalore_results_template_r1 = os.path.join(args['out_dir'], 'trim', '{sample_id}_R1.fastq.gz')
    trimgalore_results_template_r2 = os.path.join(args['out_dir'], 'trim', '{sample_id}_R2.fastq.gz')
    metrics_directory = os.path.join(args['out_dir'], 'metrics')
    fastqc_1_html_template = os.path.join(metrics_directory, 'fastqc', '{sample_id}_R1.fastqc.html')
    fastqc_1_zip_template = os.path.join(metrics_directory, 'fastqc', '{sample_id}_R1.fastqc.zip')
    fastqc_2_html_template = os.path.join(metrics_directory, 'fastqc', '{sample_id}_R2.fastqc.html')
    fastqc_2_zip_template = os.path.join(metrics_directory, 'fastqc', '{sample_id}_R2.fastqc.zip')
    metrics_summary_template = os.path.join(metrics_directory, 'summary', '{sample_id}_summary.csv')
    metrics_gcmatrix_template = os.path.join(metrics_directory, 'summary', '{sample_id}_gcmatrix.csv')
    cnmatrix_template = os.path.join(metrics_directory, 'summary', '{sample_id}_cnmatrix.csv')


    metrics_summary_filename = os.path.join(metrics_directory, 'summary', 'summary.csv')

    bam_directory = os.path.join(args['out_dir'], 'bams')
    bam_template = os.path.join(bam_directory, '{sample_id}.bam')
    bam_index_template = os.path.join(bam_directory, '{sample_id}.bam.bai')

    hmmcopy_directory = os.path.join(args['out_dir'], 'hmmcopy')
    hmmcopy_wig_template = os.path.join(hmmcopy_directory, '{sample_id}_readcount.wig')
    hmmcopy_reads_template = os.path.join(hmmcopy_directory, '{sample_id}_reads.csv')
    hmmcopy_segments_template = os.path.join(hmmcopy_directory, '{sample_id}_segments.csv')
    hmmcopy_parameters_template = os.path.join(hmmcopy_directory, '{sample_id}_parameters.csv')
    hmmcopy_posteriors_template = os.path.join(hmmcopy_directory, '{sample_id}_posteriors.csv')
    hmmcopy_hmm_metrics_template = os.path.join(hmmcopy_directory, '{sample_id}_hmm_metrics.csv')

    hmmcopy_segments_filename = os.path.join(hmmcopy_directory, 'segments.csv')
    hmmcopy_reads_filename = os.path.join(hmmcopy_directory, 'reads.csv')
    hmmcopy_hmm_metrics_filename = os.path.join(hmmcopy_directory, 'hmm_metrics.csv')



    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('fastq_1_basename', 'sample_id', axes_origin=[]),
        value=fastq_1_basenames,
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('fastq_2_basename', 'sample_id', axes_origin=[]),
        value=fastq_2_basenames,
    )


    workflow.transform(
        name='parse_sample_sheet',
        ctx={'mem': 64},
        func=single_cell_nextseq.tasks.parse_sample_sheet,
        args=(
            mgd.InputFile(args['samplesheet']),
            mgd.OutputFile(sample_info_filename),
        ),
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
        func=single_cell_nextseq.workflow.create_fastqc_workflow,
        args=(mgd.InputFile('fastq_1', 'sample_id', fnames=fastq_1_filenames),
              mgd.InputFile('fastq_2', 'sample_id', fnames=fastq_2_filenames),
              mgd.OutputFile('fastq_trim_1', 'sample_id', template=trimgalore_results_template_r1),
              mgd.OutputFile('fastq_trim_2', 'sample_id', template=trimgalore_results_template_r2),
              config,
              metrics_directory,
              mgd.InputInstance('sample_id'),
              # mgd.TempInputObj('fastq_1_basename', 'sample_id'),
              # mgd.TempInputObj('fastq_2_basename', 'sample_id'),
            ),
        )


    workflow.subworkflow(
        name='alignment_workflow',
        axes=('sample_id',),
        func=single_cell_nextseq.workflow.create_alignment_workflow,
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
        func=single_cell_nextseq.workflow.create_hmmcopy_workflow,
        args=(
            mgd.InputFile('bam', 'sample_id', template=bam_template),
            mgd.OutputFile(hmmcopy_wig_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_reads_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_segments_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_parameters_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_posteriors_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_hmm_metrics_template, 'sample_id'),
            mgd.OutputFile(cnmatrix_template, 'sample_id'),
            mgd.InputInstance('sample_id'),
            config,
        ),
    )


    workflow.subworkflow(
        name='summary_workflow',
        func=single_cell_nextseq.workflow.create_summary_workflow,
        args=(
            mgd.InputFile(hmmcopy_segments_template, 'sample_id'),
            mgd.InputFile(hmmcopy_reads_template, 'sample_id'),
            mgd.InputFile(hmmcopy_hmm_metrics_template, 'sample_id'),
            mgd.InputFile(metrics_summary_template, 'sample_id'),
            mgd.InputFile(metrics_gcmatrix_template, 'sample_id'),
            mgd.InputFile(sample_info_filename),
            mgd.InputFile(cnmatrix_template, 'sample_id'),
            config,
            args,
        ),
    )

    pyp.run(workflow)

