import yaml
import os

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


def read_samplesheet_hiseq(args):

    lib_id = os.path.normpath(args['hiseq_dir'])
    lib_id = os.path.basename(lib_id)

    fastq_1_filenames = {}
    fastq_2_filenames = {}
    
    sample_ids = []

    with open(args['samplesheet']) as freader:
        header = True
    
        for line in freader:
            line = line.strip().split(',')
    
            if header:
                if line[0] == "Description":
                    desc = line[1]
                elif line[0] == "[Data]":
                    header = False
            else:
                if line[0] in ['Sample_ID', 'Sample-ID']:
                    continue
                sampid = line[0]
                samp_idx = line[5]
                samp_idx2 = line[7]

                dirname = lib_id + '_' + samp_idx + '-' + samp_idx2
                sample_ids.append(sampid)

                for lane in args['lanes']:
                    basename = lane + '_' + samp_idx + '-' + samp_idx2
                    fq1 = os.path.join(args['hiseq_dir'], lane, dirname,
                                       basename + "_1.fastq.gz")
                    fq2 = os.path.join(args['hiseq_dir'], lane, dirname,
                                       basename + "_2.fastq.gz")
                    fastq_1_filenames[(sampid,lane)] = fq1
                    fastq_2_filenames[(sampid,lane)] = fq2                

    return desc, sample_ids, fastq_1_filenames, fastq_2_filenames

def getpath(vals):
    return os.path.join(*vals)
