import os
import argparse
import utils
import pypeliner
import pypeliner.managed as mgd
from workflows import fastq_preprocessing_workflow, alignment_workflow, hmmcopy_workflow, summary_workflow, merge_workflow, wgs_workflow, snv_postprocessing, variant_calling_workflow, realignment_workflow

from workflows import bam_postprocessing_workflow

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

    parser.add_argument('--matched_normal',
                        help='''Path to matched wgs normal.''')

    parser.add_argument('--lanes',
                        nargs='*',
                        help='''Lanes to analyze.''')

    parser.add_argument('--nextseq',
                        action='store_true',
                        help='''Lanes to analyze.''')

    parser.add_argument('--realign',
                        action='store_true',
                        help='''Lanes to analyze.''')

    parser.add_argument('--generate_pseudo_wgs',
                        action='store_true',
                        help='''Lanes to analyze.''')


    args = vars(parser.parse_args())
    args['tmpdir'] = os.path.join(args['out_dir'], 'tmp')
    
    if args['matched_normal'] and not args['generate_pseudo_wgs']:
        raise Exception('generate_pseudo_wgs must be set if matched_normal is provided')
    
    return args

def main():

    args = parse_args()

    trim = False if args['nextseq'] else True

    pyp = pypeliner.app.Pypeline(config=args)

    config = utils.load_config(args)

    desc, sample_ids, fastq_1_filenames, fastq_2_filenames = utils.read_samplesheet_hiseq(args)

    workflow = pypeliner.workflow.Workflow()

    fastq_dir = os.path.join(args['out_dir'], 'trimmed_fastq')
    trimgalore_results_template_r1 = os.path.join(fastq_dir, '{lane}', 'trim', '{sample_id}_R1.fastq.gz')
    trimgalore_results_template_r2 = os.path.join(fastq_dir, '{lane}', 'trim', '{sample_id}_R2.fastq.gz')

    bam_directory = os.path.join(args['out_dir'], 'bams')
    bam_template = os.path.join(bam_directory, '{sample_id}.bam')
    bam_index_template = os.path.join(bam_directory, '{sample_id}.bam.bai')

    pseudo_wgs_bam = os.path.join(args['out_dir'], 'pseudo_wgs', 'merged.sorted.markdups.bam')
    pseudo_wgs_bai = os.path.join(args['out_dir'], 'pseudo_wgs', 'merged.sorted.markdups.bam.bai')

    snv = os.path.join(args['out_dir'], 'pseudo_wgs', 'variant_calling', 'overlapping_calls.csv')
    countdata = os.path.join(args['out_dir'], 'pseudo_wgs', 'counts', 'counts.csv')


    vals = [(sample,lane) for lane in args['lanes'] for sample in sample_ids]
    
    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'lane'),
        value=vals,
    )


    #run fastqc and trim galore on fastq files (per lane per 
    workflow.subworkflow(
        name='fastqc_workflow',
        axes=('sample_id', 'lane',),
        func=fastq_preprocessing_workflow.create_fastq_workflow,
        args=(
              mgd.InputFile('fastq_1', 'sample_id', 'lane', fnames=fastq_1_filenames),
              mgd.InputFile('fastq_2', 'sample_id', 'lane', fnames=fastq_2_filenames),
              mgd.OutputFile('fastq_trim_1', 'sample_id', 'lane', template=trimgalore_results_template_r1),
              mgd.OutputFile('fastq_trim_2', 'sample_id', 'lane', template=trimgalore_results_template_r2),
              config,
              mgd.InputInstance('lane'),
              mgd.InputInstance('sample_id'),
              args,
              trim
            ),
        )

    workflow.subworkflow(
        name='alignment_workflow',
        axes=('sample_id', 'lane',),
        func=alignment_workflow.create_alignment_workflow,
        args=(
            mgd.InputFile('fastq_trim_1', 'sample_id', 'lane', template=trimgalore_results_template_r1),
            mgd.InputFile('fastq_trim_2', 'sample_id', 'lane', template=trimgalore_results_template_r2),
            mgd.TempOutputFile('aligned_per_cell_per_lane.sorted.bam', 'sample_id', 'lane'),
            mgd.InputFile(config['ref_genome']),
            mgd.InputInstance('lane'),
            mgd.InputInstance('sample_id'),
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
            mgd.TempInputFile('aligned_per_cell_per_lane.sorted.bam', 'sample_id', 'lane'),
            mgd.TempOutputFile('merged_lanes.bam', 'sample_id'),
            mgd.TempOutputFile('merged_lanes.bam.bai', 'sample_id'),
            config,
            args['lanes']
        ),
    )

 
    workflow.subworkflow(
        name='realignment_workflow',
        func=realignment_workflow.create_realignment_workflow,
        args=(
            mgd.TempInputFile('merged_lanes.bam', 'sample_id'),
            mgd.TempOutputFile('merged_realign.bam', 'sample_id', axes_origin=[]),
            config,
            args,
            sample_ids
        ),
    )
  
  
  
    #merge bams per sample for all lanes
    workflow.subworkflow(
        name='bam_postprocess_workflow',
        axes=('sample_id',),
        func=bam_postprocessing_workflow.create_bam_post_workflow,
        args=(
            mgd.TempInputFile('merged_realign.bam', 'sample_id'),
            mgd.OutputFile('bam_markdups', 'sample_id', template=bam_template),
            mgd.OutputFile('bam_markdups_index', 'sample_id', template=bam_index_template),
            mgd.InputFile(config['ref_genome']),
            mgd.TempOutputFile('merge_metrics_summary_template', 'sample_id'),
            mgd.TempOutputFile('merge_metrics_gcmatrix_template', 'sample_id'),
            mgd.InputInstance('sample_id'),
            mgd.InputFile(args['samplesheet']),
            config,
            args,
            desc,
            args['lanes']
        ),
    )
   
   
   
    workflow.subworkflow(
        name='hmmcopy_workflow',
        axes=('sample_id',),
        func=hmmcopy_workflow.create_hmmcopy_workflow,
        args=(
            mgd.InputFile('bam_markdups', 'sample_id', template=bam_template),
            mgd.TempOutputFile('hmmcopy_reads_template', 'sample_id'),
            mgd.TempOutputFile('hmmcopy_segments_template', 'sample_id'),
            mgd.TempOutputFile('hmmcopy_hmm_metrics_template', 'sample_id'),
            mgd.TempOutputFile('cnmatrix_template', 'sample_id'),
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
            mgd.TempInputFile('hmmcopy_segments_template', 'sample_id'),
            mgd.TempInputFile('hmmcopy_reads_template', 'sample_id'),
            mgd.TempInputFile('hmmcopy_hmm_metrics_template', 'sample_id'),
            mgd.TempInputFile('merge_metrics_summary_template', 'sample_id'),
            mgd.TempInputFile('merge_metrics_gcmatrix_template', 'sample_id'),
            mgd.TempInputFile('cnmatrix_template', 'sample_id'),
            config,
            args,
            sample_ids
        ),
    )
   
   
    if args['generate_pseudo_wgs']:
   
        workflow.subworkflow(
            name='wgs_workflow',
            func=wgs_workflow.create_wgs_workflow,
            args = (
                    mgd.InputFile('bam_markdups', 'sample_id', template=bam_template),
                    mgd.OutputFile(pseudo_wgs_bam),
                    mgd.OutputFile(pseudo_wgs_bai),
                    mgd.InputFile(config['ref_genome']),
                    sample_ids,
                    mgd.InputFile(args['samplesheet']),
                    config,
                    args,
                    desc,
            )
        )
  
  
    if args['matched_normal']:
        
        workflow.subworkflow(
                name='varcalls',
                func=variant_calling_workflow.create_varcall_workflow,
                args=(
                      mgd.InputFile(pseudo_wgs_bam),
                      mgd.InputFile(args['matched_normal']),
                      mgd.InputFile(config['ref_genome']),
                      mgd.OutputFile(snv),
                      config,
                      args
                ),
            )
        
        workflow.subworkflow(
                             name='postprocessing',
                             func=snv_postprocessing.create_snv_postprocessing_workflow,
                             args=(
                                   mgd.InputFile('bam_markdups', 'sample_id', template=bam_template),
                                   mgd.InputFile(snv),
                                   mgd.OutputFile(countdata),
                                   sample_ids,
                                   config,
                                   args
                                   )
                             )

    pyp.run(workflow)

if __name__ == '__main__':
    main()
