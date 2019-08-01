'''
Created on Apr 25, 2017

@author: dgrewal
'''

import argparse
from collections import defaultdict
import os, subprocess, math
from shutil import copyfile
import logging

class Utils(object):
    '''
    helper functions for the parsers
    in the bigdata pipeline.
    '''

    @staticmethod
    def test_args(args):
        '''
        if using the all_files option, then we dont need params, tid and nid as
        it is included in file
        '''

        if not args.infile:
            if hasattr(args, 'params_file') and args.paramsfile:
                argparse.ArgumentTypeError('Paramsfile is only required'
                                           ' when infile is provided.')

            if args.tumour_id:
                argparse.ArgumentTypeError('tumour_id is only required'
                                           ' when infile is provided.')

            if args.normal_id:
                argparse.ArgumentTypeError('normal_id is only required'
                                           ' when infile is provided.')

            if args.case_id:
                argparse.ArgumentTypeError('case_id is only required'
                                           ' when infile is provided.')


        if args.infile:
            if hasattr(args, 'params_file') and not args.paramsfile:
                argparse.ArgumentTypeError('Paramsfile is required'
                                           ' when infile is provided.')
            if args.tumour_id:
                argparse.ArgumentTypeError('tumour_id is required'
                                           ' when infile is provided.')
            if args.normal_id:
                argparse.ArgumentTypeError('normal_id is required'
                                           ' when infile is provided.')

            if args.case_id:
                argparse.ArgumentTypeError('case_id is only required'
                                           ' when infile is provided.')
                
    # pylint: disable=too-many-arguments
    @classmethod
    def get_inputs(cls, tid, nid, case, infile, all_files, paramsfile=None,
                   fh_names=None, titanfile=None):
        '''
        load the all_files if provided, o.w. use the infile and
        paramsfile params
        '''
        if not fh_names:
            fh_names = ['params', 'segments']

        output = defaultdict()
        if infile:
            if titanfile and paramsfile:
                val = (paramsfile, titanfile, infile)
            elif paramsfile:
                val = (paramsfile, infile)
            else:
                val = infile
            output[(case, tid, nid)] = val
        else:
            freader = open(all_files)

            header = freader.readline()
            if 'titan_results' in header:
                fh_names += ['titan_results']
            
            header = cls.build_indices(header)
            k_idx = [header['case_id'],
                     header['tumour_id'],
                     header['normal_id']
                     ]

            if isinstance(fh_names, str):
                v_idx = header[fh_names]
            else:
                v_idx = [header[val] for val in fh_names]

            for line in freader:
                # ignore comments
                if line[0] == '#':
                    continue
                line = line.strip().split()
                key = tuple(line[val] for val in k_idx)

                if isinstance(v_idx, list):
                    val = tuple(line[val] for val in v_idx)
                else:
                    val = line[v_idx]

                output[key] = val
        return output

    @staticmethod
    def get_label_mapping(lmap_file, prefix='label_'):
        '''
        adding label_ prefix to labels -> will help downstream to identify
        if a field is label
        '''
        if not lmap_file:
            return {}

        data = defaultdict(lambda: defaultdict(str))

        with open(lmap_file) as label_map:
            for line in label_map:
                line = line.strip().split()

                if line[0] == 'case_id':
                    labels = [prefix + val for val in line[1:]]
                    idxs = [i for i,v in enumerate(labels) if v not in 
                            [prefix+"tumour_id",prefix+"normal_id",prefix+"cohort"]]
                    labels = [labels[idx] for idx in idxs]
                    data["labels"] = labels
                    continue

                case = line[0]
                vals = line[1:]
                vals = [vals[idx] for idx in idxs]

                data[case] = vals

        return data

    @staticmethod
    def read_file_to_list(fname):
        '''
        simplistic file parser, reads file into list,
        removes line breaks
        '''

        if not fname:
            return

        data = [val.strip() for val in open(fname).readlines()]
        return data



    @staticmethod
    def build_indices(line, colnames=None):
        '''
        converts the line to dict
        with value as key and index as val
        get the indices for the provided cols
        '''

        if isinstance(line, str):
            line = line.strip().split()

        indices = dict([(val, i) for (i, val) in enumerate(line)])

        if colnames:
            col_indices = [indices[val] for val in colnames]

            return col_indices

        return indices
    
    
    @staticmethod
    def write_list(ofile, vals, labs=None, sep='\t'):
        '''
        write the list to file (join by sep)
        '''
        if labs:
            vals = vals + labs

        vals = [str(val) for val in vals]
        outstr = sep.join(vals) + '\n'
        ofile.write(outstr)

    @staticmethod
    def rm_low_mapp_regions(premapp, output, mappref, rscript):
        '''
        use Yikan's R script to remove the low mapp regions
        '''
        if not mappref:
            copyfile(premapp, output)
            return

        path_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   'filter_low_mappability_'
                                   'breakpoints_component.R')

        cmd = [rscript, path_script,
               premapp, mappref,
               output]

        cmd = ' '.join(cmd)

        proc = subprocess.Popen(cmd, shell=True)

        proc.communicate()
        retc = proc.returncode

        if retc:
            raise Exception('Unexpected error when removing low'
                            ' mappability regions')

    @staticmethod
    def get_chr_length(build='GRCh37'):
        '''
        returns lengths of all 24 chromosomes
        for a given build
        '''
        if build == 'GRCh37':
            # default chr_lengths from hg19 chromosome (build: GRCh37 released
            # February 27, 2009) please see
            # www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/index.shtml
            chr_lengths = {'1': 249250621, '2': 243199373, '3': 198022430,
                           '4': 191154276, '5': 180915260, '6': 171115067,
                           '7': 159138663, '8': 146364022, '9': 141213431,
                           '10': 135534747, '11': 135006516, '12': 133851895,
                           '13': 115169878, '14': 107349540, '15': 102531392,
                           '16': 90354753, '17': 81195210, '18': 78077248,
                           '19': 59128983, '20': 63025520, '21': 48129895,
                           '22': 51304566, 'X': 155270560, 'Y': 59373566}
        else:
            raise Exception('build %s not supported' % build)

        return chr_lengths
    
    @staticmethod
    def parse_case_order(filename):
        '''
        load the case-ordering file
        '''
        if not filename:
            return defaultdict(set)

        file_reader = open(filename)

        out = defaultdict(list)

        for line in file_reader:
            if line.strip() == '':
                continue

            line = line.strip().split('\t')
            if line[0] == 'project':
                continue

            out[line[0]] = [val.strip(' \t\n\r') for val in line[1].split(',')]

        return out
    @staticmethod
    def build_new_label_dict(filename):
        '''
        loads the label mapping file for renaming
        labels in the plots
        '''
        if not filename:
            return

        out = defaultdict(str)
        freader = open(filename)
        for line in freader:
            line = line.strip().split()
            if line[0] == 'label':
                continue
            out[line[0]] = line[1]
        return out


    @staticmethod
    def get_new_label(label, label_map_dict):
        '''
        get new label from label_dict
        '''
        if not label_map_dict:
            return label

        if isinstance(label, list):
            nlab = [label_map_dict[val] if val in label_map_dict else val
                    for val in label]
            nlab = list(set(nlab))
        else:
            if label in label_map_dict:
                nlab = label_map_dict[label]
            else:
                nlab = label
        return nlab
    
    @staticmethod
    def get_patterns_contexts():
        '''
        returns trinucelotide patterns and contexts
        '''
        contexts = defaultdict(list)
        contexts['C'] = ['ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'CCC', 'CCG',
                         'CCT', 'GCA', 'GCC', 'GCG', 'GCT', 'TCA', 'TCC',
                         'TCG', 'TCT']
        contexts['T'] = ['ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'CTC', 'CTG',
                         'CTT', 'GTA', 'GTC', 'GTG', 'GTT', 'TTA', 'TTC',
                         'TTG', 'TTT']

        sub_patterns = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']

        return sub_patterns, contexts
    
    @staticmethod
    def get_sub_pattern(ref, alt, trinuc=None):
        '''
        calculates subpattern
        '''
        def reverse_complement(base):
            '''
            returns reverse complement
            of genomic bases
            '''
            if base is None:
                return base
            out = ""
            for letter in base:
                if letter == "A":
                    out += "T"
                elif letter == "T":
                    out += "A"
                elif letter == "C":
                    out += "G"
                elif letter == "G":
                    out += "C"
                else:
                    logging.getLogger("single_cell.wgs.parse").warn('Invalid TC')
                    return
            out = out[::-1]  # reverse the string
            return out

        if len(ref) > 1 or len(alt)>1:
            return None,None

        if ref == 'A' or ref == 'G':
            ref = reverse_complement(ref)
            alt = reverse_complement(alt)

            if trinuc:
                trinuc = reverse_complement(trinuc)

        if not ref or not alt:
            return None, None

        return ref + '>' + alt, trinuc
    @staticmethod
    def get_keys(indict, keys):
        """
        goes through the keys and returns
        """
        def recurse_dict(data, outdata, i):
            '''
            recurse into dict and yield ylims.
            '''
            for key, val in data.items():
                if isinstance(val, dict):
                    outdata[i].add(key)
                    i+=1
                    i =recurse_dict(val, outdata, i)
                else:
                    outdata[i].add(key)
            return i-1
    
        i=0
        outdata = [set() for _ in range(len(keys))]
        recurse_dict(indict, outdata, i)
    
        outdata = [sorted(v) for v in outdata]
    
        return outdata
    @staticmethod
    def collapse_dict(data):
        '''
        given a dict, code tries to find the highest and lowest data points
        '''
        def recurse_dict(data, out):
            '''
            recurse into dict and yield ylims.
            '''
            for key, val in data.items():
                if isinstance(val, dict):
                    recurse_dict(val, out[key])
                else:
                    if isinstance(val, list) or isinstance(val, set):
                        out.default_factory = int
                        out[key] = len(val)

                    else:
                        raise Exception('not a list')

        getdict = lambda: defaultdict(getdict)
        out_dict = getdict()

        recurse_dict(data, out_dict)

        return out_dict
    
    # pylint: disable=R0913
    @staticmethod
    def write_data(fwriter, proj, cases, typ, vals, typ_suffix=''):
        '''
        writes the plot data diles
        '''
        if isinstance(cases, list) and isinstance(vals, list):
            for case, val in zip(cases, vals):
                ostr = '\t'.join(
                    [proj, case, typ + typ_suffix, str(val)]) + '\n'
                fwriter.write(ostr)

        elif isinstance(vals, list) and isinstance(typ, list):
            for tval, val in zip(typ, vals):
                try:
                    val = '0' if math.isnan(val) else val
                except:
                    pass
                ostr = '\t'.join(
                    [proj, cases, tval + typ_suffix, str(val)]) + '\n'
                fwriter.write(ostr)

        elif isinstance(vals, list):
            for val in vals:
                ostr = '\t'.join(
                    [proj, cases, typ + typ_suffix, str(val)]) + '\n'
                fwriter.write(ostr)

        elif isinstance(typ, list):
            for tval in typ:
                ostr = '\t'.join(
                    [proj, cases, tval + typ_suffix, str(vals)]) + '\n'
                fwriter.write(ostr)

        else:
            ostr = '\t'.join([proj, cases, typ + typ_suffix, str(vals)]) + '\n'
            fwriter.write(ostr)
