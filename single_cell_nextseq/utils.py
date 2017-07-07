import yaml
import os
import warnings

def load_config(args):
    try:
        with open(args['config_file']) as infile:
            config = yaml.load(infile)

    except IOError:
        print 'Unable to open config file: {0}'.format(args['config_file'])
    return config

def read_samplesheet_nextseq(args, fastq_directory):
    try:
        with open(args['samplesheet']) as infile:
            lines = [x.strip('\n').strip(',') for x in infile.readlines()]
        
        run_id = [s.split(',')[1] for s in lines if 'Experiment Name,' in s][0]
        
        library_id = [s.split(',')[1] for s in lines if 'Description,' in s][0]
        
        start_index = lines.index('[Data]')+2

        num_samples = len(lines[start_index:])

        sample_ids = []
        fastq_1_filenames = {}
        fastq_2_filenames = {}
        for i, line in zip(range(num_samples), lines[start_index:]):
            sample_id = line.split(',')[0]
            sample_id = sample_id.replace('-', '_')
            r1_basename = '{0}_S{1}_R1_001'.format(sample_id, str(i+1))
            r2_basename = '{0}_S{1}_R2_001'.format(sample_id, str(i+1))
            fastq_1_filenames[sample_id] = os.path.join(fastq_directory, r1_basename + '.fastq.gz')
            fastq_2_filenames[sample_id] = os.path.join(fastq_directory, r2_basename + '.fastq.gz')
            sample_ids.append(sample_id)

    except IOError:
        print 'Unable to open file \'SampleSheet.csv\' in directory: {0}'.format(args['nextseq_dir'])

    return run_id, library_id, sample_ids, fastq_1_filenames, fastq_2_filenames



def read_samplesheet_hiseq(args, fastq_directory):
    try:
        with open(args['samplesheet']) as infile:
            lines = [x.strip('\n').strip(',') for x in infile.readlines()]
        
        run_id = [s.split(',')[1] for s in lines if 'Experiment Name,' in s][0]
        
        library_id = [s.split(',')[1] for s in lines if 'Description,' in s][0]
        
        start_index = lines.index('[Data]')+2

        num_samples = len(lines[start_index:])

        sample_ids = []
        fastq_1_filenames = {}
        fastq_2_filenames = {}
        for i, line in zip(range(num_samples), lines[start_index:]):
            sample_id = line.split(',')[0]
            sample_id = sample_id.replace('-', '_')
            r1_basename = '{0}_S{1}_R1_001'.format(sample_id, str(i+1))
            r2_basename = '{0}_S{1}_R2_001'.format(sample_id, str(i+1))
            fastq_1_filenames[sample_id] = os.path.join(fastq_directory, r1_basename + '.fastq.gz')
            fastq_2_filenames[sample_id] = os.path.join(fastq_directory, r2_basename + '.fastq.gz')
            sample_ids.append(sample_id)

    except IOError:
        print 'Unable to open file \'SampleSheet.csv\' in directory: {0}'.format(args['nextseq_dir'])

    return run_id, library_id, sample_ids, fastq_1_filenames, fastq_2_filenames


def get_readgroup_template(library_id, run_id, config):
    if 'read_group' in config.keys():
        read_group_template = (
            '@RG\tID:' + str(library_id) + '_{sample_id}_' + str(run_id) + 
            '\tPL:' + str(config['read_group']['PL']) +
            '\tPU:' + str(run_id) +
            '\tLB:' + str(library_id) + '_{sample_id}' +
            '\tSM:' + '{sample_id}' +
            '\tCN:' + str(config['read_group']['CN']))
    
    else:
        warnings.warn('Config file does not contain read group information! ' + 
                      'This will affect duplicate marking if BAMs are later merged. ' +
                      'Creating BAM without read group information in header.')

    return read_group_template

def getpath(vals):
    return os.path.join(*vals)
