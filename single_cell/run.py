import os
import argparse
import utils
import itertools
import pypeliner
import pypeliner.managed as mgd
from workflows import hmmcopy
from workflows import strelka
from workflows import alignment
from workflows import merge_bams
from workflows import pseudo_wgs
from workflows import realignment
from workflows import mutationseq
from workflows import singlecell_summary
from workflows import snv_postprocessing
from workflows import alignment_postprocessing

import pypeliner.storage


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    pypeliner.app.add_arguments(parser)

    parser.add_argument('sample_info',
                        help='''Per sample meta data CSV''')

    parser.add_argument('fastqs_file',
                        help='''Path to input fastq table CSV.''')

    parser.add_argument('library_id',
                        help='''Library id.''')

    parser.add_argument('out_dir',
                        help='''Path to output files.''')
 
    parser.add_argument('config_file',
                        help='''Path to yaml config file.''')

    parser.add_argument('--matched_normal',
                        help='''Path to matched wgs normal.''')

    parser.add_argument('--realign',
                        action='store_true',
                        help='''Lanes to analyze.''')

    parser.add_argument('--generate_pseudo_wgs',
                        action='store_true',
                        help='''Lanes to analyze.''')

    args = vars(parser.parse_args())
    args['tmpdir'] = os.path.join(args['out_dir'], 'tmp')

    if args['matched_normal'] and not args['generate_pseudo_wgs']:
        raise Exception('generate_pseudo_wgs must be'
                        ' set if matched_normal is provided')

    return args


def main():
    
    storage = pypeliner.storage.AzureBlobStorage()

    args = parse_args()

    pyp = pypeliner.app.Pypeline(config=args)

    config = utils.load_config(args)

    fastq1_files, fastq2_files, sampleids, lanes, seqinfo = utils.read_fastqs_file(args)

    for fileset in (fastq1_files, fastq2_files):
        for (cell_id, lane_id), filename in fileset.iteritems():
            blob_name = filename.replace('/datadrive/analysis/', '')
            assert 'datadrive' not in blob_name
            store = storage.create_store(blob_name)
            if not store.get_exists():
                print 'uploading', blob_name
                store.allocated_filename = filename
                store.push()
            else:
                print 'exists', blob_name
            fileset[(cell_id, lane_id)] = blob_name

    blob_name = args['sample_info'].replace('/datadrive/analysis/', '')
    store = storage.create_store(blob_name)
    if not store.get_exists():
        store.allocated_filename = args['sample_info']
        store.push()
    args['sample_info'] = blob_name

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'lane'),
        value=fastq1_files.keys(),
    )

    workflow.subworkflow(
        name='alignment_workflow',
        axes=('sample_id', 'lane',),
        func=alignment.create_alignment_workflow,
        args=(
            mgd.InputFile('fastq_1', 'sample_id', 'lane', fnames=fastq1_files),
            mgd.InputFile('fastq_2', 'sample_id', 'lane', fnames=fastq2_files),
            mgd.TempOutputFile('aligned_per_cell_per_lane.sorted.bam', 'sample_id', 'lane'),
            config['ref_genome'],
            mgd.InputInstance('lane'),
            mgd.InputInstance('sample_id'),
            config,
            args,
            seqinfo
        ),
    )

    # merge bams per sample for all lanes
    workflow.subworkflow(
        name='merge_workflow',
        axes=('sample_id',),
        func=merge_bams.create_merge_workflow,
        args=(
            mgd.TempInputFile('aligned_per_cell_per_lane.sorted.bam', 'sample_id', 'lane'),
            mgd.TempOutputFile('merged_lanes.bam', 'sample_id'),
            mgd.TempOutputFile('merged_lanes.bam.bai', 'sample_id'),
            config,
            lanes,
        ),
    )

    workflow.subworkflow(
        name='realignment_workflow',
        func=realignment.create_realignment_workflow,
        args=(
            mgd.TempInputFile('merged_lanes.bam', 'sample_id'),
            mgd.TempInputFile('merged_lanes.bam.bai', 'sample_id'),
            mgd.TempOutputFile('merged_realign.bam', 'sample_id', axes_origin=[]),
            config,
            args['out_dir'],
            args['realign'],
            sampleids
        ),
    )

    bam_directory = os.path.join(args['out_dir'], 'bams')
    bam_template = os.path.join(bam_directory, '{sample_id}.bam')
    bam_index_template = os.path.join(bam_directory, '{sample_id}.bam.bai')
    results_dir = os.path.join(args['out_dir'], 'results')
    gc_metrics = os.path.join(results_dir, '{}_gc_metrics.csv'.format(args['library_id']))
    # merge bams per sample for all lanes
    workflow.subworkflow(
        name='bam_postprocess_workflow',
        func=alignment_postprocessing.create_bam_post_workflow,
        args=(
            mgd.TempInputFile('merged_realign.bam', 'sample_id'),
            mgd.OutputFile('bam_markdups', 'sample_id', template=bam_template, axes_origin=[]),
            mgd.OutputFile('bam_markdups_index', 'sample_id',
                           template=bam_index_template, axes_origin=[]),
            config['ref_genome'],
            mgd.TempOutputFile('alignment_metrics.csv'),
            mgd.OutputFile(gc_metrics),
            sampleids,
            config,
            args['out_dir'],
        ),
    )

    results_dir = os.path.join(args['out_dir'], 'results')
    segs_filename = os.path.join(results_dir, '{}_segments.csv'.format(args['library_id']))
    reads_filename = os.path.join(results_dir, '{}_reads.csv'.format(args['library_id']))
    workflow.subworkflow(
        name='hmmcopy_workflow',
        func=hmmcopy.create_hmmcopy_workflow,
        args=(
            mgd.InputFile('bam_markdups', 'sample_id', template=bam_template),
            mgd.InputFile('bam_markdups_index', 'sample_id', template=bam_index_template),
            mgd.OutputFile(reads_filename),
            mgd.OutputFile(segs_filename),
            mgd.TempOutputFile('hmmcopy_hmm_metrics.csv'),
            sampleids,
            config,
            args
        ),
    )
  
    # merge all samples per lane together
    workflow.subworkflow(
        name='summary_workflow',
        func=singlecell_summary.create_summary_workflow,
        args=(
            mgd.InputFile(args['sample_info']),
            mgd.InputFile(segs_filename),
            mgd.InputFile(reads_filename),
            mgd.TempInputFile('hmmcopy_hmm_metrics.csv'),
            mgd.TempInputFile('alignment_metrics.csv'),
            mgd.InputFile(gc_metrics),
            config,
            args,
            sampleids
        ),
    )
  
    if args['generate_pseudo_wgs']:
  
        pseudo_wgs_bam = os.path.join(args['out_dir'], 'pseudo_wgs',
                                      'merged.sorted.markdups.bam')
        pseudo_wgs_bai = os.path.join(args['out_dir'], 'pseudo_wgs',
                                      'merged.sorted.markdups.bam.bai')
        workflow.subworkflow(
            name='wgs_workflow',
            func=pseudo_wgs.create_wgs_workflow,
            args=(
                mgd.InputFile('bam_markdups', 'sample_id',
                              template=bam_template),
                mgd.OutputFile(pseudo_wgs_bam),
                mgd.OutputFile(pseudo_wgs_bai),
                config['ref_genome'],
                sampleids,
                config,
                args['out_dir'],
            )
        )
  
    if args['matched_normal']:
        varcalls_dir = os.path.join(args['out_dir'], 'pseudo_wgs',
                                    'variant_calling')
        museq_vcf = os.path.join(varcalls_dir, 'museq_snv.vcf')
        museq_csv = os.path.join(varcalls_dir, 'museq_snv.csv')
        workflow.subworkflow(
            name='museq',
            func=mutationseq.create_museq_workflow,
            args=(
                mgd.InputFile(pseudo_wgs_bam),
                mgd.InputFile(args['matched_normal']),
                config['ref_genome'],
                mgd.OutputFile(museq_vcf),
                mgd.OutputFile(museq_csv),
                config,
                args,
            ),
        )
  
        strelka_snv_vcf = os.path.join(varcalls_dir, 'strelka_snv.vcf')
        strelka_indel_vcf = os.path.join(varcalls_dir, 'strelka_indel.vcf')
        strelka_snv_csv = os.path.join(varcalls_dir, 'strelka_snv.csv')
        strelka_indel_csv = os.path.join(varcalls_dir, 'strelka_indel.csv')
        workflow.subworkflow(
            name='strelka',
            func=strelka.create_strelka_workflow,
            args=(
                mgd.InputFile(pseudo_wgs_bam),
                mgd.InputFile(args['matched_normal']),
                config['ref_genome'],
                mgd.OutputFile(strelka_indel_vcf),
                mgd.OutputFile(strelka_snv_vcf),
                mgd.OutputFile(strelka_indel_csv),
                mgd.OutputFile(strelka_snv_csv),
            ),
        )
  
        countdata = os.path.join(args['out_dir'], 'pseudo_wgs',
                                 'counts', 'counts.csv')
        workflow.subworkflow(
            name='postprocessing',
            func=snv_postprocessing.create_snv_postprocessing_workflow,
            args=(
                mgd.InputFile('bam_markdups', 'sample_id',
                              template=bam_template),
                mgd.InputFile(museq_csv),
                mgd.InputFile(strelka_snv_csv),
                mgd.OutputFile(countdata),
                sampleids,
                config,
                args['out_dir'],
            )
        )
 
    pyp.run(workflow)

if __name__ == '__main__':
    main()
