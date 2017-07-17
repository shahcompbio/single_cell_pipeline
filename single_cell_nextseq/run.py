import os
import argparse
import utils
import pypeliner
import pypeliner.managed as mgd
from workflows import fastqc_workflow, alignment_workflow, hmmcopy_workflow, summary_workflow, merge_workflow

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    pypeliner.app.add_arguments(parser)

    parser.add_argument('hiseq_dir',
                        help='''Path to input hiseq directory.''')

    parser.add_argument('samplesheet',
                        help='''Path to the samplesheet.''')

    parser.add_argument('out_dir',
                        help='''Path to output files.''')

    parser.add_argument('config_file',
                        help='''Path to yaml config file.''')

    parser.add_argument('--lanes',
                        nargs='*',
                        help='''Lanes to analyze.''')

    args = vars(parser.parse_args())
    args['tmpdir'] = os.path.join(args['out_dir'], 'tmp')

    return args


def main():

    args = parse_args()

    pyp = pypeliner.app.Pypeline(config=args)

    config = utils.load_config(args)

    lanes, desc, sample_ids, fastq_1_filenames, fastq_2_filenames = utils.read_samplesheet_hiseq(args)

    workflow = pypeliner.workflow.Workflow()


    perlane_dir = os.path.join(args['out_dir'], 'lanes')

    trimgalore_results_template_r1 = os.path.join(perlane_dir, '{lane}', 'trim', '{sample_id}_R1.fastq.gz')
    trimgalore_results_template_r2 = os.path.join(perlane_dir, '{lane}', 'trim', '{sample_id}_R2.fastq.gz')

    bam_template = os.path.join(perlane_dir, '{lane}', 'bams', '{sample_id}.bam')
    bam_index_template = os.path.join(perlane_dir, '{lane}', 'bams', '{sample_id}.bam.bai')
    metrics_summary_template = os.path.join(perlane_dir, '{lane}', 'summary', '{sample_id}_summary.csv')
    metrics_gcmatrix_template = os.path.join(perlane_dir, '{lane}', 'summary', '{sample_id}_gcmatrix.csv')

    bam_directory = os.path.join(args['out_dir'], 'bams')
    merge_template = os.path.join(bam_directory, '{sample_id}.merged.bam')
    merge_index_template = os.path.join(bam_directory, '{sample_id}.merged.bam.bai')

    metrics_directory = os.path.join(args['out_dir'], 'metrics')
    merge_metrics_summary_template = os.path.join(metrics_directory, 'summary', '{sample_id}_summary_all_lanes.csv')
    merge_metrics_gcmatrix_template = os.path.join(metrics_directory, 'summary', '{sample_id}_gcmatrix_all_lanes.csv')

    hmmcopy_directory = os.path.join(args['out_dir'], 'hmmcopy')
    hmmcopy_reads_template = os.path.join(hmmcopy_directory, '{sample_id}_reads.csv')
    hmmcopy_segments_template = os.path.join(hmmcopy_directory, '{sample_id}_segments.csv')
    hmmcopy_hmm_metrics_template = os.path.join(hmmcopy_directory, '{sample_id}_hmm_metrics.csv')
    cnmatrix_template = os.path.join(hmmcopy_directory, '{sample_id}_cnmatrix.csv')


    vals = [(sample,lane) for lane in lanes for sample in sample_ids]
    
    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'lane'),
        value=vals,
    )


    workflow.subworkflow(
        name='fastqc_workflow',
        axes=('sample_id', 'lane',),
        func=fastqc_workflow.create_fastqc_workflow,
        args=(
              mgd.InputFile('fastq_1', 'sample_id', 'lane', fnames=fastq_1_filenames),
              mgd.InputFile('fastq_2', 'sample_id', 'lane', fnames=fastq_2_filenames),
              mgd.OutputFile('fastq_trim_1', 'sample_id', 'lane', template=trimgalore_results_template_r1),
              mgd.OutputFile('fastq_trim_2', 'sample_id', 'lane', template=trimgalore_results_template_r2),
              config,
              mgd.InputInstance('lane'),
              mgd.InputInstance('sample_id'),
              args,
              True
            ),
        )

    workflow.subworkflow(
        name='alignment_workflow',
        axes=('sample_id', 'lane',),
        func=alignment_workflow.create_alignment_workflow,
        args=(
            mgd.InputFile('fastq_trim_1', 'sample_id', 'lane', template=trimgalore_results_template_r1),
            mgd.InputFile('fastq_trim_2', 'sample_id', 'lane', template=trimgalore_results_template_r2),
            mgd.OutputFile('bam', 'sample_id', 'lane', template=bam_template),
            mgd.OutputFile('bam_index', 'sample_id', 'lane', template=bam_index_template),
            mgd.InputFile(config['ref_genome']),
            mgd.OutputFile(metrics_summary_template, 'sample_id','lane'),
            mgd.OutputFile(metrics_gcmatrix_template, 'sample_id','lane'),
            mgd.InputInstance('lane'),
            mgd.InputInstance('sample_id'),
            mgd.InputFile(args['samplesheet']),
            config,
            args,
            desc
        ),
    )

    #merge bams per sample for all lanes
    workflow.subworkflow(
        name='merge_workflow',
        axes=('sample_id',),
        func=merge_workflow.create_merge_workflow,
        args=(
            mgd.InputFile('bam', 'sample_id', 'lane', template=bam_template),
            mgd.OutputFile('bam', 'sample_id', template=merge_template),
            mgd.OutputFile('bam_index', 'sample_id', template=merge_index_template),
            mgd.InputFile(config['ref_genome']),
            mgd.OutputFile(merge_metrics_summary_template, 'sample_id'),
            mgd.OutputFile(merge_metrics_gcmatrix_template, 'sample_id'),
            mgd.InputInstance('sample_id'),
            mgd.InputFile(args['samplesheet']),
            config,
            args,
            desc,
            lanes
        ),
    )
  
    workflow.subworkflow(
        name='hmmcopy_workflow',
        axes=('sample_id',),
        func=hmmcopy_workflow.create_hmmcopy_workflow,
        args=(
            mgd.InputFile('bam', 'sample_id', template=merge_template),
            mgd.OutputFile(hmmcopy_reads_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_segments_template, 'sample_id'),
            mgd.OutputFile(hmmcopy_hmm_metrics_template, 'sample_id'),
            mgd.OutputFile(cnmatrix_template, 'sample_id'),
            mgd.InputInstance('sample_id'),
            config,
            args
        ),
    )

    # merge all samples per lane together
    workflow.subworkflow(
        name='summary_workflow',
        func=summary_workflow.create_summary_workflow,
        args=(
            mgd.InputFile(hmmcopy_segments_template, 'sample_id'),
            mgd.InputFile(hmmcopy_reads_template, 'sample_id'),
            mgd.InputFile(hmmcopy_hmm_metrics_template, 'sample_id'),
            mgd.InputFile(merge_metrics_summary_template, 'sample_id'),
            mgd.InputFile(merge_metrics_gcmatrix_template, 'sample_id'),
            mgd.InputFile(cnmatrix_template, 'sample_id'),
            config,
            args,
            sample_ids
        ),
    )
 
    pyp.run(workflow)

if __name__ == '__main__':
    main()
