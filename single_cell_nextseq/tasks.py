import pypeliner.commmandline


def demultiplex_fastq_files(sample_sheet_filename, fastq_filenames, nextseq_directory, temp_directory):
    out_dir = os.path.join(config['out_dir'], 'fastq')
    
    make_parent_directory(out_dir)

    pypeliner.commmandline(
        'bcl2fastq',
        '--runfolder-dir=' + nextseq_directory, 
        '--output-dir=' + temp_directory, 
        '--no-lane-splitting')


    cmd = config['bcl2fastq']
    
    cmd_args = [
                '--runfolder-dir=' + nextseq_dir, 
                '--output-dir=' + out_dir, 
                '--no-lane-splitting'
                ]
    
    run_cmd(cmd, cmd_args, max_mem=34, queue=config['queue'], hosts=config['hosts'])

