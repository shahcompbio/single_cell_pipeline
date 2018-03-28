'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd

import tasks


def create_copyclone_workflow(bam_file, bai_file, reads, segments, metrics, sample_ids, config, args, results_dir ):

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )

    workflow.transform(
        name='correct_reads',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.correct_reads,
        axes=('sample_id',),
        args=(
            mgd.InputFile('bam_markdups', 'sample_id', fnames=bam_file),
            mgd.InputFile('bai_markdups', 'sample_id', fnames=bai_file),
            mgd.TempOutputFile('corrected_reads.csv', 'sample_id'),
            config,
            mgd.TempSpace('hmmcopy_temp', 'sample_id'),
            mgd.InputInstance("sample_id")
        ),
    )

    workflow.transform(
        name='merge_corrected_reads',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.merge_reads,
        args=(
            mgd.TempInputFile('corrected_reads.csv', 'sample_id'),
            mgd.TempOutputFile('corrected_reads.csv'),
        ),
    )

    workflow.transform(
        name="run_copyclone",
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.run_copyclone,
        args=(
            mgd.TempInputFile('corrected_reads.csv'),
            mgd.OutputFile(reads),
            mgd.OutputFile(segments),
            mgd.OutputFile(metrics),
        ),
        kwargs={
            "A": config["copyclone"]["A"],
            "alpha_A": config["copyclone"]["alpha_A"],
            "pi": config["copyclone"]["pi"],
            "alpha_pi": config["copyclone"]["alpha_pi"],
            "tau": config["copyclone"]["tau"],
            "nu": config["copyclone"]["nu"],
            "eta": config["copyclone"]["eta"],
            "shape": config["copyclone"]["shape"],
            "rate": config["copyclone"]["rate"],
            "ploidy_states": config["copyclone"]["ploidy_states"],
            "num_states": config["copyclone"]["num_states"],
        }
    )

    return workflow
