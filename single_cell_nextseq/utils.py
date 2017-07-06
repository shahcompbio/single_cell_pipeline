import yaml
import os


def load_config(args):
    try:
        with open(args['config_file']) as file:
            config = yaml.load(file)

    except IOError as e:
        print 'Unable to open config file: {0}'.format(args['config_file'])
    return config

def read_samplesheet(args, fastq_directory):
    try:
        with open(args['samplesheet']) as file:
            lines = [x.strip('\n').strip(',') for x in file.readlines()]
        
        run_id = [s.split(',')[1] for s in lines if 'Experiment Name,' in s][0]
        
        library_id = [s.split(',')[1] for s in lines if 'Description,' in s][0]
        
        start_index = lines.index('[Data]')+2

        num_samples = len(lines[start_index:])

        sample_ids = []
        fastq_1_basenames = {}
        fastq_2_basenames = {}
        fastq_1_filenames = {}
        fastq_2_filenames = {}
        for i, line in zip(range(num_samples), lines[start_index:]):
            sample_id = line.split(',')[0]
            sample_id = sample_id.replace('-', '_')
            fastq_1_basenames[sample_id] = '{0}_S{1}_R1_001'.format(sample_id, str(i+1))
            fastq_2_basenames[sample_id] = '{0}_S{1}_R2_001'.format(sample_id, str(i+1))
            fastq_1_filenames[sample_id] = os.path.join(fastq_directory, fastq_1_basenames[sample_id] + '.fastq.gz')
            fastq_2_filenames[sample_id] = os.path.join(fastq_directory, fastq_2_basenames[sample_id] + '.fastq.gz')
            sample_ids.append(sample_id)

    except IOError as e:
        print 'Unable to open file \'SampleSheet.csv\' in directory: {0}'.format(args['nextseq_dir'])

    return run_id, library_id, sample_ids, fastq_1_filenames, fastq_2_filenames, fastq_1_basenames, fastq_2_basenames

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
