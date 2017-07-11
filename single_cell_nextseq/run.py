import os
import argparse
import utils
import pypeliner
import pypeliner.managed as mgd
from workflows import tasks, fastqc_workflow, alignment_workflow, hmmcopy_workflow, summary_workflow



scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
plot_hmmcopy_script = os.path.join(scripts_directory, 'plot_hmmcopy.py')



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

    args = vars(parser.parse_args())
    args['tmpdir'] = os.path.join(args['out_dir'], 'tmp')

    return args


def merge_bams(inputs, output):
    print inputs
    open(output, 'w').close()


def main():

    args = parse_args()

    pyp = pypeliner.app.Pypeline(config=args)

    config = utils.load_config(args)

    lanes, desc, sample_ids = utils.read_samplesheet_hiseq(args)


    fastq_r1 = os.path.join(args['out_dir'], 'fastq', '{lane}','{sample_id}_R1.fastq.gz')
    fastq_r2 = os.path.join(args['out_dir'], 'fastq', '{lane}','{sample_id}_R2.fastq.gz')

    trimgalore_results_template_r1 = os.path.join(args['out_dir'], 'trim', '{lane}','{sample_id}_R1.fastq.gz')
    trimgalore_results_template_r2 = os.path.join(args['out_dir'], 'trim', '{lane}','{sample_id}_R2.fastq.gz')

    bam_directory = os.path.join(args['out_dir'], 'bams')
    bam_template = os.path.join(bam_directory, '{lane}', '{sample_id}.bam')
    bam_index_template = os.path.join(bam_directory, '{lane}', '{sample_id}.bam.bai')
    merge_template = os.path.join(bam_directory, '{sample_id}.merged.bam')
    merge_index_template = os.path.join(bam_directory, '{sample_id}.merged.bam.bai')



    metrics_directory = os.path.join(args['out_dir'], 'metrics')
    metrics_summary_template = os.path.join(metrics_directory, 'summary', '{lane}', '{sample_id}_summary.csv')
    metrics_gcmatrix_template = os.path.join(metrics_directory, 'summary', '{lane}', '{sample_id}_gcmatrix.csv')

    merge_metrics_summary_template = os.path.join(metrics_directory, 'summary', '{sample_id}_summary_all_lanes.csv')
    merge_metrics_gcmatrix_template = os.path.join(metrics_directory, 'summary', '{sample_id}_gcmatrix_all_lanes.csv')



    hmmcopy_directory = os.path.join(args['out_dir'], 'hmmcopy')
    hmmcopy_reads_template = os.path.join(hmmcopy_directory, '{lane}', '{sample_id}_reads.csv')
    hmmcopy_segments_template = os.path.join(hmmcopy_directory, '{lane}', '{sample_id}_segments.csv')
    hmmcopy_hmm_metrics_template = os.path.join(hmmcopy_directory, '{lane}', '{sample_id}_hmm_metrics.csv')
    cnmatrix_template = os.path.join(metrics_directory, 'summary', '{lane}', '{sample_id}_cnmatrix.csv')



    workflow = pypeliner.workflow.Workflow()

    vals = [(lane,sample) for lane in lanes for sample in sample_ids]
    
    
    workflow.setobj(
        obj=mgd.OutputChunks('lane', 'sample_id'),
        value=vals,
    )

    workflow.transform(
                       name='get_fastq_files',
                       axes=('lane', 'sample_id',),
                       func=tasks.get_fastq_files,
                       args=(mgd.InputFile(args['samplesheet']),
                             mgd.InputFile(args['hiseq_dir']),
                             mgd.OutputFile('fastq_1', 'lane', 'sample_id', template=fastq_r1),
                             mgd.OutputFile('fastq_2', 'lane', 'sample_id', template=fastq_r2),
                             mgd.InputInstance('sample_id'),
                             mgd.InputInstance('lane'),
                             )
                       )


    workflow.subworkflow(
        name='fastqc_workflow',
        axes=('lane', 'sample_id',),
        func=fastqc_workflow.create_fastqc_workflow,
        args=(
              mgd.InputFile('fastq_1', 'lane', 'sample_id', template=fastq_r1),
              mgd.InputFile('fastq_2', 'lane', 'sample_id', template=fastq_r2),
              mgd.OutputFile('fastq_trim_1', 'lane', 'sample_id', template=trimgalore_results_template_r1),
              mgd.OutputFile('fastq_trim_2', 'lane', 'sample_id', template=trimgalore_results_template_r2),
              config,
              mgd.InputInstance('lane'),
              mgd.InputInstance('sample_id'),
              args,
              False
            ),
        )

    workflow.subworkflow(
        name='alignment_workflow',
        axes=('lane', 'sample_id',),
        func=alignment_workflow.create_alignment_workflow,
        args=(
            mgd.InputFile('fastq_trim_1', 'lane', 'sample_id', template=trimgalore_results_template_r1),
            mgd.InputFile('fastq_trim_2', 'lane', 'sample_id', template=trimgalore_results_template_r2),
            mgd.OutputFile('bam', 'lane', 'sample_id', template=bam_template),
            mgd.OutputFile('bam_index', 'lane', 'sample_id', template=bam_index_template),
            mgd.InputFile(config['ref_genome']),
            mgd.OutputFile(metrics_summary_template, 'lane','sample_id'),
            mgd.OutputFile(metrics_gcmatrix_template, 'lane','sample_id'),
            mgd.InputInstance('lane'),
            mgd.InputInstance('sample_id'),
            mgd.InputFile(args['samplesheet']),
            config,
            args,
            desc
        ),
    )

    workflow.subworkflow(
        name='hmmcopy_workflow',
        axes=('lane', 'sample_id',),
        func=hmmcopy_workflow.create_hmmcopy_workflow,
        args=(
            mgd.InputFile('bam', 'lane', 'sample_id', template=bam_template),
            mgd.OutputFile(hmmcopy_reads_template, 'lane', 'sample_id'),
            mgd.OutputFile(hmmcopy_segments_template, 'lane', 'sample_id'),
            mgd.OutputFile(hmmcopy_hmm_metrics_template, 'lane', 'sample_id'),
            mgd.OutputFile(cnmatrix_template, 'lane', 'sample_id'),
            mgd.InputInstance('lane'),
            mgd.InputInstance('sample_id'),
            config,
            args
        ),
    )


    # merge all samples per lane together
    workflow.subworkflow(
        name='summary_workflow',
        func=summary_workflow.create_summary_workflow,
        axes=('lane',),
        args=(
            mgd.InputFile(hmmcopy_segments_template, 'lane', 'sample_id'),
            mgd.InputFile(hmmcopy_reads_template, 'lane', 'sample_id'),
            mgd.InputFile(hmmcopy_hmm_metrics_template, 'lane', 'sample_id'),
            mgd.InputFile(metrics_summary_template, 'lane', 'sample_id'),
            mgd.InputFile(metrics_gcmatrix_template, 'lane', 'sample_id'),
            mgd.InputFile(cnmatrix_template, 'lane', 'sample_id'),
            config,
            args,
            mgd.InputInstance('lane'),
            sample_ids
        ),
    )

    
    #testing merge by samples
    workflow.transform(
        name='merge_bams',
        axes = ('sample_id',),
        func=merge_bams,
        args=(
            mgd.InputFile(hmmcopy_segments_template, 'lane', 'sample_id'),
            mgd.TempOutputFile('merged.bam', 'sample_id'),
        ),
    )


    #merge bams per sample for all lanes
#     workflow.subworkflow(
#         name='merge_workflow',
#         axes=('sample_id',),
#         func=merge_workflow.create_merge_workflow,
#         args=(
#             mgd.InputFile('bam', 'lane', 'sample_id', template=bam_template),
#             mgd.OutputFile('bam', 'sample_id', template=merge_template),
#             mgd.OutputFile('bam_index', 'sample_id', template=merge_index_template),
#             mgd.InputFile(config['ref_genome']),
#             mgd.OutputFile(merge_metrics_summary_template, 'sample_id'),
#             mgd.OutputFile(merge_metrics_gcmatrix_template, 'sample_id'),
#             'all',
#             mgd.InputInstance('sample_id'),
#             mgd.InputFile(args['samplesheet']),
#             config,
#             args,
#             desc
#         ),
#     )




    pyp.run(workflow)

if __name__ == '__main__':
    main()

