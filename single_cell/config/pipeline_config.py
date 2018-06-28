'''
Created on Jun 6, 2018

@author: dgrewal
'''

import single_cell
from config_reference import extract_from_reference
import yaml
import copy


def get_version():
    version = single_cell.__version__
    # strip setuptools metadata
    version = version.split("+")[0]
    return version


def get_config_params(override=None):
    version = get_version()
    input_params = {
        "cluster": "azure", "aligner": "bwa-aln",
        "reference": "grch37", "smoothing_function": "modal",
        "bin_size": 500000, "copynumber_bin_size": 1000,
        "version": version
    }

    if override:
        input_params.update(override)

    input_params["version"] = input_params["version"].replace('.', '_')

    return input_params


def get_pools(reference, version):

    pools = {"standard": "singlecell{}standard_{}".format(reference, version),
             "highmem": "singlecell{}highmem_{}".format(reference, version),
             "multicore": "singlecell{}multicore_{}".format(reference, version)
             }

    return {"pools": pools}


def get_picard_wgs_params():

    picard_params = {"min_bqual": extract_from_reference(["min_bqual"]),
                     "min_mqual": extract_from_reference(["min_mqual"]),
                     "count_unpaired": extract_from_reference(['count_unpaired'])
                     }
    return {'picard_wgs_params': picard_params}


def get_titan_params(cluster, reference, binsize):

    chrom_info_filename = extract_from_reference(
        [cluster, reference, 'chrom_info_filename'])

    window_size = binsize

    ref_data_dir = extract_from_reference(
        [cluster, reference, 'copynumber_ref_data'])

    # Binned GC content file
    gc_wig = extract_from_reference(
        [cluster, reference, 'gc_wig_file', binsize])

    # Binned mappability file
    mappability_wig = extract_from_reference(
        [cluster, reference, 'map_wig_file', binsize])

    titan_params = {"normal_contamination": extract_from_reference(['normal_contamination']),
                    'num_clusters': extract_from_reference(['num_clusters']),
                    'ploidy': extract_from_reference(['ploidy']),
                    'chrom_info_filename': chrom_info_filename,
                    'window_size': window_size,
                    'ref_data_dir': ref_data_dir,
                    'gc_wig': gc_wig,
                    'mappability_wig': mappability_wig
                    }

    return {'titan_params': titan_params}


def get_hmmcopy_params(cluster, reference, binsize, smoothing_function):

    params = {
        'multipliers': extract_from_reference(['multipliers']),
        'map_cutoff': extract_from_reference(['map_cutoff']),
        'bin_size': binsize,
        'e': extract_from_reference(['e']),
        'eta': extract_from_reference(['eta']),
        'g': extract_from_reference(['g']),
        'lambda': extract_from_reference(['lambda']),
        'min_mqual': extract_from_reference(['min_mqual']),
        'nu': extract_from_reference(['nu']),
        'num_states': extract_from_reference(['num_states_hmmcopy']),
        's': extract_from_reference(['s']),
        'strength': extract_from_reference(['strength']),
        'kappa': extract_from_reference(['kappa']),
        'm': extract_from_reference(['m']),
        'mu': extract_from_reference(['mu']),
        'smoothing_function': smoothing_function,
        'exclude_list': extract_from_reference([cluster, reference, 'exclude_list']),
        'gc_wig_file': extract_from_reference([cluster, reference, 'gc_wig_file', binsize]),
        'map_wig_file': extract_from_reference([cluster, reference, 'map_wig_file', binsize]),
    }

    return {"hmmcopy_params": {"autoploidy": params}}

def get_copyclone_params(cluster, reference, binsize, smoothing_function):
    ploidy_2_params = {'A':[0.994, 0.994, 0.9994, 0.994, 0.994, 0.994, 0.994],
                       'alpha_A': [1000, 1000, 10000, 1000, 1000, 1000, 1000],
                       'pi':[0.05, 0.1, 0.5, 0.2, 0.05, 0.05, 0.05],
                       'alpha_pi':[2, 2, 50, 2, 2, 2, 2]}

    ploidy_3_params = {'A':[0.994, 0.994, 0.994, 0.9994, 0.994, 0.994, 0.994],
                       'alpha_A': [1000, 1000, 1000, 10000, 1000, 1000, 1000],
                       'pi':[0.05, 0.05, 0.1, 0.5, 0.2, 0.05, 0.05],
                       'alpha_pi':[2, 2, 2, 50, 2, 2, 2]}


    ploidy_4_params = {'A':[0.994, 0.994, 0.994, 0.994, 0.9994, 0.994, 0.994],
                       'alpha_A': [1000, 1000, 1000, 1000, 10000, 1000, 1000],
                       'pi':[0.05, 0.05, 0.05, 0.1, 0.5, 0.2, 0.05],
                       'alpha_pi':[2, 2, 2, 2, 50, 2, 2]}



    ploidy_states = {2: ploidy_2_params, 3: ploidy_3_params, 4:ploidy_4_params}

    params = {
        'map_cutoff': extract_from_reference(['map_cutoff']),
        'bin_size': binsize,
        'gc_wig_file': extract_from_reference([cluster, reference, 'gc_wig_file', binsize]),
        'map_wig_file': extract_from_reference([cluster, reference, 'map_wig_file', binsize]),
        'smoothing_function': smoothing_function,
        'exclude_list': None,
        'min_mqual': extract_from_reference(['min_mqual']),
        'num_states': extract_from_reference(['num_states_copyclone']),
        'tau': [500, 25, 25, 25, 25, 25, 15],
        'nu': [5, 5, 5, 5, 5, 5, 5],
        'eta': [5000, 5000, 5000, 5000, 5000, 5000, 5000],
        'shape': [3, 30, 30, 30, 30, 30, 20],
        'rate': [0.01, 1, 1, 1, 1, 1, 1],
        'ploidy_states': ploidy_states,
    }

    return {"copyclone": params}


def get_cell_filter():

    return {"good_cells": extract_from_reference(['good_cells'])}


def get_memory_requests():

    params = {'memory': extract_from_reference(['memory'])}

    return params


def get_global_params(cluster, reference, aligner):

    params = {
        # Reference genome to align to in fasta format. Should be indexed by
        # aligner and samtools
        'ref_genome': extract_from_reference([cluster, reference, 'ref_genome']),
        # Chromosomes for copy number analysis
        'chromosomes': extract_from_reference([cluster, reference, 'chromosomes']),
        'mutationseq_python': extract_from_reference([cluster, 'mutationseq_python']),
        'mutationseq': extract_from_reference([cluster, 'mutationseq']),
        'mutationseq_model': extract_from_reference([cluster, 'mutationseq_model']),
        'gc_windows': extract_from_reference([cluster, reference, 'gc_windows']),
        'one_split_job': extract_from_reference([cluster, 'one_split_job']),
        'max_cores': extract_from_reference([cluster, 'max_cores']),
        'aligner': aligner,
        'adapter': extract_from_reference(['adapter']),
        'adapter2': extract_from_reference(['adapter2']),
        'split_size': extract_from_reference(['split_size']),
    }

    return params


def get_destruct_params(cluster, reference):

    params = {
        'destruct': {
            'genome_fasta': extract_from_reference([cluster, reference, 'ref_genome']),
            'genome_fai': extract_from_reference([cluster, reference, 'ref_genome_index']),
        }
    }

    return params


def get_databases():

    status_data = {
        'kwargs': {
            'split_size': extract_from_reference(['split_size'])
        }
    }

    databases = {
        "databases": extract_from_reference(['databases']),
        'cosmic_status': copy.deepcopy(status_data),
        'dbsnp_status': copy.deepcopy(status_data),
        'mappability': copy.deepcopy(status_data),
        'snpeff_status': copy.deepcopy(status_data),
        'tri_nucleotide_context': copy.deepcopy(status_data),
    }

    return databases


def get_singlecell_pipeline_config(config_params):

    reference = config_params["reference"]
    cluster = config_params["cluster"]

    params = {}

    params.update(get_pools(reference, config_params["version"]))

    params.update(get_picard_wgs_params())

    params.update(
        get_titan_params(
            cluster, reference, config_params["copynumber_bin_size"]
        )
    )

    params.update(
        get_hmmcopy_params(
            cluster, reference, config_params["bin_size"],
            config_params["smoothing_function"]
        )
    )

    params.update(
        get_copyclone_params(
            cluster, reference, config_params["bin_size"],
            config_params["smoothing_function"]
        )
    )

    params.update(
        get_global_params(
            cluster, reference, config_params["aligner"]
        )
    )

    params.update(get_cell_filter())

    params.update(get_memory_requests())

    params.update(get_destruct_params(cluster, reference))

    params.update(get_databases())

    return params


def write_config(params, filepath):
    with open(filepath, 'w') as outputfile:
        yaml.safe_dump(params, outputfile, default_flow_style=False)
