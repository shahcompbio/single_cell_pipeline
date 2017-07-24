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
                                       museq_parsed,
                                       strelka_parsed,
                                       output,
                                       sample_ids,
                                       config,
                                       args
                                       ):
    script_path = os.path.join(os.path.realpath(os.path.dirname(__file__)),
                               'scripts', 'get_cell_counts.py')

    countdata = os.path.join(args['out_dir'], 'pseudo_wgs', 'counts', '{sample_id}_counts.csv')

    scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
    merge_tables_script = os.path.join(scripts_directory, 'merge.py')



    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )



    workflow.commandline(
                       name='overlap_var_calls',
                       ctx={'mem': config['med_mem']},
                       args=(
                            config['python'],
                            merge_tables_script,
                            '--merge_type', 'inner',
                            '--nan_value', 'NA',
                            '--input', mgd.InputFile(museq_parsed), mgd.InputFile(strelka_parsed),
                            '--key_cols', 'case_id', 'chromosome', 'start', 'stop', 'ref', 'alt',
                            '--separator', 'tab',
                            '--type', 'merge',
                            '--output', mgd.TempOutputFile("overlapping_calls.csv"),
                             )
                       )



    workflow.commandline(
        name='count_reads',
        axes=('sample_id',),
        ctx={'mem': config['low_mem']},
        args=(
              config['python'],
              script_path,
              mgd.InputFile('bam', 'sample_id', fnames=bam_file),
              mgd.TempInputFile("overlapping_calls.csv"),
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