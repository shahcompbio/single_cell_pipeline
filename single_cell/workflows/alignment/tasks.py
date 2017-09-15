import pypeliner.commandline
import os

def get_readgroup(run_id, sample_id, args, config, seqinfo):
    platform = 'illumina'
    centre = 'UBCBRC' if seqinfo[sample_id] == 'nextseq' else 'BCCAGSC'

    if 'read_group' in config:
        if config['read_group']['PL']:
            platform = str(config['read_group']['PL'])
        if config['read_group']['CN']:
            centre = str(config['read_group']['CN'])
    
    library_id = args['library_id']
    read_group_template = (
        '@RG\tID:' + str(library_id) + '_' + sample_id + '_' + str(run_id) +
        '\tPL:' + platform +
        '\tPU:' + str(run_id) +
        '\tLB:' + str(library_id) + '_' + sample_id +
        '\tSM:' + sample_id +
        '\tCN:' + centre)

    return read_group_template


def bam_sort(bam_filename, sorted_bam_filename, config):
    pypeliner.commandline.execute(
        'picard', '-Xmx12G',
        'SortSam',
        'INPUT=' + bam_filename,
        'OUTPUT=' + sorted_bam_filename,
        'SORT_ORDER=coordinate',
        'VALIDATION_STRINGENCY=LENIENT',
        'MAX_RECORDS_IN_RAM=5000000')


def align_paired_end(fastq1, fastq2, output, tempdir,
                     reference, config, readgroup):
    """
    run bwa aln on both fastq files,
    bwa sampe to align, and convert to bam with samtools view
    """
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    read_1_sai = os.path.join(tempdir, 'read_1.sai')
    read_2_sai = os.path.join(tempdir, 'read_2.sai')



    pypeliner.commandline.execute(
            'bwa',
            'aln',
            reference,
            fastq1,
            '>',
            read_1_sai                                  
      )



    pypeliner.commandline.execute(
            'bwa',
            'aln',
            reference,
            fastq2,
            '>',
            read_2_sai,
                                  
      )

    pypeliner.commandline.execute(
            'bwa', 'sampe',
            '-r', readgroup,
            reference,
            read_1_sai,
            read_2_sai,
            fastq1,
            fastq2,
            '|',
            'samtools', 'view',
            '-bSh', '-',
            '>',
            output,
        )

    
