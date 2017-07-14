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
                                       snv_calls,
                                       output,
                                       sample_ids,
                                       config,
                                       args
                                       ):
    script_path = os.path.join(os.path.realpath(os.path.dirname(__file__)),
                               'scripts', 'get_cell_counts.py')

    countdata = os.path.join(args['out_dir'], 'counts', '{sample_id}_counts.csv')

    
    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )


    workflow.commandline(
        name='count_reads',
        axes=('sample_id',),
        ctx={'mem': config['low_mem']},
        args=(
              config['python'],
              script_path,
              mgd.InputFile('bam', 'sample_id', fnames=bam_file),
              mgd.InputFile(snv_calls),
              mgd.OutputFile(countdata, 'sample_id'),
              mgd.InputInstance('sample_id')
        ),
    )

    workflow.transform(
        name='merge_counts',
        ctx={'mem': config['low_mem']},
        func=tasks.merge_csv,
        args=(
              mgd.InputFile(countdata, 'sample_id'),
              mgd.OutputFile(output),
              'outer',
              'chrom,coord,ref_base,var_base'
        ),
    )



    return workflow