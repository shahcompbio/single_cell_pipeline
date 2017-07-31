import pypeliner.commandline
import warnings
import os

def get_readgroup(library_id, run_id, config, sample_id):
    if 'read_group' in config.keys():
        read_group_template = (
            '@RG\tID:' + str(library_id) + '_' + sample_id + '_' + str(run_id) +
            '\tPL:' + str(config['read_group']['PL']) +
            '\tPU:' + str(run_id) +
            '\tLB:' + str(library_id) + '_' + sample_id +
            '\tSM:' + sample_id +
            '\tCN:' + str(config['read_group']['CN']))

    else:
        warnings.warn('Config file does not contain read group information! ' +
                      'This will affect duplicate marking if BAMs are later merged. ' +
                      'Creating BAM without read group information in header.')

    return read_group_template


def bam_sort(bam_filename, sorted_bam_filename, config):
    pypeliner.commandline.execute(
        config['picard'], '-Xmx12G',
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

    read_1_sai = os.path.join(tempdir, 'read_1.sai')
    read_2_sai = os.path.join(tempdir, 'read_2.sai')



    pypeliner.commandline.execute(
            config['bwa'],
            'aln',
            reference,
            fastq1,
            '>',
            read_1_sai,
                                  
      )

    pypeliner.commandline.execute(
            config['bwa'],
            'aln',
            reference,
            fastq2,
            '>',
            read_2_sai,
                                  
      )

    pypeliner.commandline.execute(
            config['bwa'], 'sampe',
            '-r', readgroup,
            reference,
            read_1_sai,
            read_2_sai,
            fastq1,
            fastq2,
            '|',
            config['samtools'], 'view',
            '-bSh', '-',
            '>',
            output,
        )

    
