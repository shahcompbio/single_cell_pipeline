'''
Created on Jul 14, 2017

@author: dgrewal


take snv calls, bams from cells and get counts for ref and alt
'''
import os
import tasks
import pypeliner
import pypeliner.managed as mgd


def create_snv_postprocessing_workflow(
    bam_file,
    bai_file,
    museq_parsed,
    strelka_parsed,
    output,
    overlapping_calls,
    sample_ids,
    config,
):

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )

    workflow.transform(
        name='overlap_var_calls',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.merge_csv,
        args=([mgd.InputFile(museq_parsed),
               mgd.InputFile(strelka_parsed)],
              mgd.OutputFile(overlapping_calls),
              'outer',
              ['case_id', 'chromosome', 'start', 'stop', "ref", "alt"],
        ),
        kwargs={'sep': ',', 'suffixes': ['_museq', "_strelka"]},
    )

    workflow.transform(
        name='count_reads',
        axes=('sample_id',),
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.get_counts,
        args=(
            mgd.InputFile('bam', 'sample_id', fnames=bam_file),
            mgd.InputFile('bai', 'sample_id', fnames=bai_file),
            mgd.InputFile(overlapping_calls),
            mgd.TempOutputFile("counts.csv", 'sample_id'),
            mgd.InputInstance('sample_id')
        ),
    )

    workflow.transform(
        name='merge_counts',
        ctx={'mem': config["memory"]['low'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.concat_csv,
        args=(
            mgd.TempInputFile("counts.csv", 'sample_id'),
            mgd.OutputFile(output),
        ),
    )

    return workflow
