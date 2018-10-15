'''
Created on Jun 6, 2018

@author: dgrewal
'''
import collections
import single_cell
from config_reference import extract_from_reference
import yaml
import copy
import os

def get_version():
    version = single_cell.__version__
    # strip setuptools metadata
    version = version.split("+")[0]
    return version


def get_config_params(override=None):
    version = get_version()
    input_params = {
        "cluster": "azure", "aligner": "bwa-mem",
        "reference": "grch37", "smoothing_function": "modal",
        "bin_size": 500000, "copynumber_bin_size": 1000,
        "version": version,
        'memory': {'high': 18, 'med': 6, 'low': 2}
    }

    input_params = override_config(input_params, override)

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
                    'mappability_wig': mappability_wig,
                    'chromosomes': map(str, range(1,23)) + ['X'],
                    }

    return {'titan_params': titan_params}


def get_hmmcopy_params(cluster, reference, binsize, smoothing_function):

    params = {
        'multipliers': [1, 2, 3, 4, 5, 6],
        'map_cutoff': extract_from_reference(['map_cutoff']),
        'bin_size': binsize,
        'e': 0.999999,
        'eta': 50000,
        'g': 3,
        'lambda': 20,
        'min_mqual': extract_from_reference(['min_mqual']),
        'nu': 2.1,
        'num_states': 12,
        's': 1,
        'strength': 1000,
        'kappa': '100,100,700,100,25,25,25,25,25,25,25,25',
        'm': '0,1,2,3,4,5,6,7,8,9,10,11',
        'mu': '0,1,2,3,4,5,6,7,8,9,10,11',
        'smoothing_function': smoothing_function,
        'exclude_list': extract_from_reference([cluster, reference, 'exclude_list']),
        'gc_wig_file': extract_from_reference([cluster, reference, 'gc_wig_file', binsize]),
        'map_wig_file': extract_from_reference([cluster, reference, 'map_wig_file', binsize]),
        'classifier_training_data': extract_from_reference([cluster, 'classifier_training_data'])
    }

    return {"hmmcopy_params": {"autoploidy": params}}


def get_copyclone_params(cluster, reference, binsize, smoothing_function):
    params = {
        'map_cutoff': extract_from_reference(['map_cutoff']),
        'bin_size': binsize,
        'gc_wig_file': extract_from_reference([cluster, reference, 'gc_wig_file', binsize]),
        'map_wig_file': extract_from_reference([cluster, reference, 'map_wig_file', binsize]),
        'smoothing_function': smoothing_function,
        'exclude_list': None,
        'min_mqual': extract_from_reference(['min_mqual']),
        'num_states': 7,
        'A': [0.994, 0.994, 0.994, 0.994, 0.994, 0.994, 0.994],
        'alpha_A': [1000, 1000, 1000, 1000, 1000, 1000, 1000],
        'alpha_pi': [2, 2, 50, 2, 2, 2, 2],
        'pi': [0.05, 0.1, 0.5, 0.2, 0.05, 0.05, 0.05],
        'tau': [500, 25, 25, 25, 25, 25, 15],
        'nu': [5, 5, 5, 5, 5, 5, 5],
        'eta': [5000, 5000, 5000, 5000, 5000, 5000, 5000],
        'shape': [3, 30, 30, 30, 30, 30, 20],
        'rate': [0.01, 1, 1, 1, 1, 1, 1],
        'ploidy_states': [2, 3, 4],
    }

    return {"copyclone": params}


def get_cell_filter():

    return {"good_cells": extract_from_reference(['good_cells'])}


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
        'snpeff': copy.deepcopy(status_data),
        'tri_nucleotide_context': copy.deepcopy(status_data),
    }

    return databases


def get_container_images(cluster):
    images = {
        'bwa': 'bwa:v0.0.1', 'samtools': 'samtools:v0.0.1',
        'python_base': 'python_base:v0.0.1', 'picard': 'picard:v0.0.1',
        'single_cell_pipeline': 'single_cell_pipeline:v0.0.2',
        'gatk': 'gatk:v0.0.1', 'fastqc': 'fastqc:v0.0.1',
        'hmmcopy': 'hmmcopy:v0.0.1', 'aneufinder': 'aneufinder:v0.0.1',
        'strelka': 'strelka:v0.0.1', 'mutationseq': 'mutationseq:v0.0.1',
        'vcftools': 'vcftools:v0.0.1', 'snpeff': 'vcftools:v0.0.1',
        'titan': 'titan:v0.0.1'}

    if cluster == 'azure':
        username = os.environ["CLIENT_ID"]
        password = os.environ["SECRET_KEY"]
        server = "shahlab.azurecr.io"

        image_urls = {k:"shahlab.azurecr.io/scp/{}:v0.0.1".format(v)
                      for k, v in images.iteritems()}
    else:
        username = None
        password = None
        server = None
        image_path = '/refdata/single_cell_pipeline.simg'
        image_urls = {v: image_path for v in images}

    image_data = {}
    for name, url in image_urls.iteritems():
        image_data[name] = {'image': url, 'server': server,
                            'username': username, 'password': password}
    images = {"images": image_data}

    return images


def get_container_params(cluster,):
    params = {}
    params.update(get_container_images(cluster))

    container_type = 'docker' if cluster == 'azure' else 'singularity'
    params["container_type"] = container_type
    params["mounts"] = ['/refdata', '/datadrive']
    if cluster == 'azure':
        params["mounts"].append('/mnt')
    return {"containers": params}


def override_config(config, override):
    def update(d, u):
        for k, v in u.iteritems():
            if isinstance(v, collections.Mapping):
                d[k] = update(d.get(k, {}), v)
            else:
                d[k] = v
        return d

    if not override:
        return config

    cfg = update(config, override)

    return cfg


def get_singlecell_pipeline_config(config_params, override=None):

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
        get_global_params(
            cluster, reference, config_params["aligner"]
        )
    )

    params.update(get_cell_filter())

    params.update({'memory': config_params['memory']})

    params.update(get_destruct_params(cluster, reference))

    params.update(get_databases())

    params.update(get_container_params(cluster))

    params = override_config(params, override)

    return params


def write_config(params, filepath):
    with open(filepath, 'w') as outputfile:
        yaml.safe_dump(params, outputfile, default_flow_style=False)
