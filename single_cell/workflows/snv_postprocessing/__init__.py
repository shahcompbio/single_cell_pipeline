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
    variant_csv_files,
    output,
    overlapping_calls,
    cell_ids,
    config,
):

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cell_ids,
    )

    workflow.transform(
        name='overlap_var_calls',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.merge_csv,
        args=([mgd.InputFile(csv_filename) for csv_filename in variant_csv_files],
              mgd.OutputFile(overlapping_calls),
              'outer',
              ['case_id', 'chromosome', 'start', 'stop', "ref", "alt"],
        ),
        kwargs={'sep': ',', 'suffixes': ['_museq', "_strelka"]},
    )

    workflow.transform(
        name='count_reads',
        axes=('cell_id',),
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.get_counts,
        args=(
            mgd.InputFile('bam', 'cell_id', fnames=bam_file),
            mgd.InputFile('bai', 'cell_id', fnames=bai_file),
            mgd.InputFile(overlapping_calls),
            mgd.TempOutputFile("counts.csv", 'cell_id'),
            mgd.InputInstance('cell_id')
        ),
    )

    workflow.transform(
        name='merge_counts',
        ctx={'mem': config["memory"]['low'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.concat_csv,
        args=(
            mgd.TempInputFile("counts.csv", 'cell_id'),
            mgd.OutputFile(output),
        ),
    )

    return workflow
