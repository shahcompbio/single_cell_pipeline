'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner
import warnings

from scripts import ParseMuseq


def run_museq(tumour, normal, reference, museq_dir, out, log, interval, config):

    script = os.path.join(museq_dir, 'classify.py')
    model = os.path.join(museq_dir, 'model_v4.1.2_anaconda_sk_0.13.1.npz')

    conf = os.path.join(museq_dir, 'metadata.config')

    cmd = ['python', script, 'normal:'+ normal, 'tumour:'+ tumour,
           'reference:'+ reference, 'model:'+ model, '--out', out,
           '--log', log, '--config', conf]

    if interval:
        if "_" in interval:
            interval = interval.split('_')
            interval = "{}:{}-{}".format(*interval)
        cmd.extend(['--interval', interval])

    pypeliner.commandline.execute(*cmd)

def concatenate_vcf(infiles, outfile):
    def get_header(infile):
        '''
        extract header from the file
        '''
        header = []
        for line in infile:
            if line.startswith('##'):
                header.append(line)
            elif line.startswith('#'):
                header.append(line)
                return header
            else:
                raise Exception('invalid header: missing #CHROM line')

        warnings.warn("One of the input files is empty")
        return []

    with open(outfile, 'w') as ofile:
        header = None

        for _,ifile in infiles.iteritems():

            with open(ifile) as f:

                if not header:
                    header = get_header(f)

                    for line in header:
                        ofile.write(line)
                else:
                    if not get_header(f) == header:
                        warnings.warn(
                            'merging vcf files with mismatching headers')

                for l in f:
                    print >> ofile, l,

def parse_museq(infile, output):
    parser = ParseMuseq(infile=infile, tid='NA', nid='NA', output=output,
                        keep_dbsnp=True,keep_1000gen=True,
                        remove_duplicates=True)
    
    parser.main()
