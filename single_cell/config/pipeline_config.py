'''
Created on Jun 6, 2018

@author: dgrewal
'''
import collections
import copy

import yaml
from single_cell.config import config_reference


def override_config(config, override):
    def update(d, u):
        for k, v in u.items():
            if isinstance(v, collections.Mapping):
                d[k] = update(d.get(k, {}), v)
            else:
                d[k] = v
        return d

    if not override:
        return config

    cfg = update(config, override)

    return cfg


def get_config_params(override=None):
    input_params = {
        "cluster": "azure", "aligner": "bwa-mem", "refdir": None,
        "reference": "grch37", "smoothing_function": "modal",
        "bin_size": 500000, "copynumber_bin_size": 1000,
        'memory': {'high': 16, 'med': 6, 'low': 2},
        'version': None
    }

    input_params = override_config(input_params, override)

    return input_params


def write_config(params, filepath):
    with open(filepath, 'w') as outputfile:
        yaml.safe_dump(params, outputfile, default_flow_style=False)


def get_hmmcopy_params(reference_dir, reference, binsize, smoothing_function, version):
    referencedata = config_reference.get_cluster_reference_data(reference_dir, reference)

    docker_containers = config_reference.containers(version)['docker']
    docker_containers = {
        'single_cell_pipeline': docker_containers['single_cell_pipeline'],
        'hmmcopy': docker_containers['hmmcopy']
    }

    params = {
        'multipliers': [1, 2, 3, 4, 5, 6],
        'map_cutoff': 0.9,
        'bin_size': binsize,
        'e': 0.999999,
        'eta': 50000,
        'g': 3,
        'lambda': 20,
        'min_mqual': 20,
        'nu': 2.1,
        'num_states': 12,
        's': 1,
        'strength': 1000,
        'kappa': '100,100,700,100,25,25,25,25,25,25,25,25',
        'm': '0,1,2,3,4,5,6,7,8,9,10,11',
        'mu': '0,1,2,3,4,5,6,7,8,9,10,11',
        'smoothing_function': smoothing_function,
        'exclude_list': referencedata['exclude_list'],
        'gc_wig_file': referencedata['gc_wig_file'][binsize],
        'map_wig_file': referencedata['map_wig_file'][binsize],
        'chromosomes': referencedata['chromosomes'],
        'ref_genome': referencedata['ref_genome'],
        'docker': docker_containers,
        'igv_segs_quality_threshold': 0.75,
        'memory': {'med': 6},
        'good_cells': [
            ['median_hmmcopy_reads_per_bin', 'ge', 50],
            ['is_contaminated', 'in', ['False', 'false', False]],
        ]
    }

    return {"hmmcopy": params}


def get_align_params(reference_dir, reference, aligner, version):
    referencedata = config_reference.get_cluster_reference_data(reference_dir, reference)
    refdata_callback = config_reference.get_cluster_reference_data

    docker_containers = config_reference.containers(version)['docker']
    docker_containers = {
        'single_cell_pipeline': docker_containers['single_cell_pipeline'],
        'fastqc': docker_containers['fastqc'],
        'samtools': docker_containers['samtools'],
        'bwa': docker_containers['bwa'],
        'picard': docker_containers['picard'],
        'trimgalore': docker_containers['trimgalore'],
        'fastq_screen': docker_containers['fastq_screen'],
    }

    params = {
        'ref_genome': referencedata['ref_genome'],
        'docker': docker_containers,
        'memory': {'med': 6},
        'aligner': aligner,
        'adapter': 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC',
        'adapter2': 'CTGTCTCTTATACACATCTGACGCTGCCGACGA',
        'picard_wgs_params': {
            "min_bqual": 20,
            "min_mqual": 20,
            "count_unpaired": False,
        },
        'chromosomes': referencedata['chromosomes'],
        'gc_windows': referencedata['gc_windows'],
        'fastq_screen_params': {
            'strict_validation': True,
            'filter_contaminated_reads': False,
            'aligner': 'bwa',
            'genomes': [
                {'name': 'grch37', 'path': refdata_callback(reference_dir, 'grch37')['ref_genome']},
                {'name': 'mm10', 'path': refdata_callback(reference_dir, 'mm10')['ref_genome']},
                {'name': 'salmon', 'path': refdata_callback(reference_dir, 'GCF_002021735')['ref_genome']},
            ]
        }
    }

    return {"alignment": params}


def get_annotation_params(reference_dir, reference, version):
    referencedata = config_reference.get_cluster_reference_data(reference_dir, reference)

    docker_containers = config_reference.containers(version)['docker']
    docker_containers = {
        'single_cell_pipeline': docker_containers['single_cell_pipeline'],
        'cell_cycle_classifier': docker_containers['cell_cycle_classifier'],
        'corrupt_tree': docker_containers['corrupt_tree']
    }

    params = {
        'docker': docker_containers,
        'memory': {'med': 6},
        'classifier_training_data': referencedata['classifier_training_data'],
        'reference_gc': referencedata['reference_gc_qc'],
        'chromosomes': referencedata['chromosomes'],
        'num_states': 12,
        'map_cutoff': 0.9,
        'ref_type': reference,
        'corrupt_tree_params': {
            'neighborhood_size': 2,
            'lower_fraction': 0.05,
            'engine_nchains': 1,
            'engine_nscans': 10000,
            'model_fpr_bound': 0.1,
            'model_fnr_bound': 0.5
        },
        'good_cells': [
            ['quality', 'ge', 0.75],
            ['experimental_condition', 'notin', ["NTC", "NCC", "gDNA", "GM"]],
            ['cell_call', 'in', ['C1']],
            ['is_contaminated', 'in', ['False', 'false', False]],
        ]
    }

    return {"annotation": params}


def get_aneufinder_params(reference_dir, reference, version):
    referencedata = config_reference.get_cluster_reference_data(reference_dir, reference)

    docker_containers = config_reference.containers(version)['docker']
    params = {
        'memory': {'med': 6},
        'docker': {
            'single_cell_pipeline': docker_containers['single_cell_pipeline'],
            'aneufinder': docker_containers['aneufinder'],
        },
        'chromosomes': referencedata['chromosomes'],
        'ref_genome': referencedata['ref_genome']
    }

    return {'aneufinder': params}


def get_merge_bams_params(reference_dir, reference, cluster, version):
    referencedata = config_reference.get_cluster_reference_data(reference_dir, reference)

    if cluster in ['azure', 'aws']:
        one_split_job = True
    else:
        one_split_job = False

    docker_containers = config_reference.containers(version)['docker']
    params = {
        'memory': {'low': 4, 'med': 6, 'high': 16},
        'max_cores': 8,
        'docker': {
            'single_cell_pipeline': docker_containers['single_cell_pipeline'],
            'samtools': docker_containers['samtools']
        },
        'ref_genome': referencedata['ref_genome'],
        'split_size': 10000000,
        'chromosomes': referencedata['chromosomes'],
        'one_split_job': one_split_job
    }
    return {'merge_bams': params}


def get_split_bam_params(reference_dir, reference, version):
    referencedata = config_reference.get_cluster_reference_data(reference_dir, reference)

    docker_containers = config_reference.containers(version)['docker']

    params = {
        'memory': {'low': 4, 'med': 6, 'high': 16},
        'max_cores': 8,
        'docker': {
            'single_cell_pipeline': docker_containers['single_cell_pipeline'],
            'samtools': docker_containers['samtools']
        },
        'ref_genome': referencedata['ref_genome'],
        'split_size': 10000000,
        'chromosomes': referencedata['chromosomes'],
        'one_split_job': True
    }

    return {'split_bam': params}


def get_germline_calling_params(reference_dir, reference, version):
    referencedata = config_reference.get_cluster_reference_data(reference_dir, reference)

    docker_containers = config_reference.containers(version)['docker']

    params = {
        'memory': {'low': 4, 'med': 6, 'high': 16},
        'max_cores': 8,
        'docker': {
            'single_cell_pipeline': docker_containers['single_cell_pipeline'],
            'samtools': docker_containers['samtools'],
            'vcftools': docker_containers['vcftools'],
            'snpeff': docker_containers['snpeff'],
        },
        'ref_genome': referencedata['ref_genome'],
        'chromosomes': referencedata['chromosomes'],
        'split_size': 10000000,
        'databases': {
            'mappability': {
                'url': 'http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/release3'
                       '/wgEncodeCrgMapabilityAlign50mer.bigWig',
                'local_path': referencedata['databases']['mappability']['local_path'],
            },
            'snpeff': {
                "db": 'GRCh37.75',
                "data_dir": referencedata['databases']['snpeff']['local_path']
            },
        },
    }

    return {'germline_calling': params}


def get_variant_calling_params(reference_dir, reference, version):
    referencedata = config_reference.get_cluster_reference_data(reference_dir, reference)

    docker_containers = config_reference.containers(version)['docker']

    status_data = {
        'kwargs': {
            'split_size': 10000000
        }
    }

    params = {
        'memory': {'low': 4, 'med': 6, 'high': 16},
        'max_cores': 8,
        'docker': {
            'single_cell_pipeline': docker_containers['single_cell_pipeline'],
            'vcftools': docker_containers['vcftools'],
            'strelka': docker_containers['strelka'],
            'mutationseq': docker_containers['mutationseq'],
        },
        'ref_genome': referencedata['ref_genome'],
        'chromosomes': referencedata['chromosomes'],
        'use_depth_thresholds': True,
        'split_size': 10000000,
        'cosmic_status': copy.deepcopy(status_data),
        'dbsnp_status': copy.deepcopy(status_data),
        'mappability': copy.deepcopy(status_data),
        'snpeff': copy.deepcopy(status_data),
        'tri_nucleotide_context': copy.deepcopy(status_data),
        'databases': {
            'cosmic': {
                'download_method': 'sftp',
                'user_name': 'awm3@sfu.ca',
                'password': 'shahlabith',
                'host': 'sftp-cancer.sanger.ac.uk',
                'remote_paths': {
                    'coding': '/files/grch37/cosmic/v75/VCF/CosmicCodingMuts.vcf.gz',
                    'non_coding': '/files/grch37/cosmic/v75/VCF/CosmicNonCodingVariants.vcf.gz',
                },
                'local_path': referencedata['databases']['cosmic']['local_path'],
            },
            'dbsnp': {
                'url': 'ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b146_GRCh37p13/VCF/common_all_20151104.vcf.gz',
                'local_path': referencedata['databases']['dbsnp']['local_path'],
            },
            'mappability': {
                'url': 'http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/release3'
                       '/wgEncodeCrgMapabilityAlign50mer.bigWig',
                'local_path': referencedata['databases']['mappability']['local_path'],
            },
            'ref_genome': {
                'url': 'http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa',
                'local_path': referencedata['ref_genome'],
            },
            'snpeff': {
                "db": 'GRCh37.75',
                "data_dir": referencedata['databases']['snpeff']['local_path']
            },
        },
        'museq_params': {
            'threshold': 0.5,
            'verbose': True,
            'purity': 70,
            'coverage': 4,
            'buffer_size': '2G',
            'mapq_threshold': 10,
            'indl_threshold': 0.05,
            'normal_variant': 25,
            'tumour_variant': 2,
            'baseq_threshold': 10,
        }

    }

    return {'variant_calling': params}


def get_copy_number_calling_params(reference_dir, reference, binsize, version):
    referencedata = config_reference.get_cluster_reference_data(reference_dir, reference)

    docker_containers = config_reference.containers(version)['docker']

    params = {
        'memory': {'low': 4, 'med': 6, 'high': 16},
        'ref_genome': referencedata['ref_genome'],
        'chromosomes': referencedata['chromosomes'],
        'split_size': 10000000,
        'max_cores': None,
        'chromosomes': referencedata['chromosomes'],
        'extract_seqdata': {},
        'ref_data_dir': referencedata['copynumber_ref_data'],
        'docker': {
            'single_cell_pipeline': docker_containers['remixt'],
            'titan': docker_containers['titan']
        },
        'titan_params': {
            "normal_contamination": [0.2, 0.4, 0.6, 0.8],
            'num_clusters': [1, 2],
            'ploidy': [1, 2, 3, 4],
            'chrom_info_filename': referencedata['chrom_info_filename'],
            'window_size': binsize,
            'gc_wig': referencedata['gc_wig_file'][binsize],
            'mappability_wig': referencedata['gc_wig_file'][binsize],
        }
    }

    return {'copy_number_calling': params}


def get_infer_haps_params(reference_dir, reference, version):
    referencedata = config_reference.get_cluster_reference_data(reference_dir, reference)

    docker_containers = config_reference.containers(version)['docker']

    params = {
        'memory': {'low': 4, 'med': 6, 'high': 16},
        'max_cores': None,
        'chromosomes': referencedata['chromosomes'],
        'extract_seqdata': {
            'genome_fasta_template': referencedata['ref_genome'],
            'genome_fai_template': referencedata['ref_genome'] + '.fai',
        },
        'ref_data_dir': referencedata['copynumber_ref_data'],
        'docker': {
            'single_cell_pipeline': docker_containers['single_cell_pipeline'],
            'remixt': docker_containers['remixt'],
        },
    }

    return {'infer_haps': params}


def get_count_haps_params(reference_dir, reference, version):
    referencedata = config_reference.get_cluster_reference_data(reference_dir, reference)

    docker_containers = config_reference.containers(version)['docker']

    params = {
        'memory': {'low': 4, 'med': 6, 'high': 16},
        'max_cores': None,
        'chromosomes': referencedata['chromosomes'],
        'extract_seqdata': {
            'genome_fasta_template': referencedata['ref_genome'],
            'genome_fai_template': referencedata['ref_genome'] + '.fai',
        },
        'ref_data_dir': referencedata['copynumber_ref_data'],
        'docker': {
            'single_cell_pipeline': docker_containers['single_cell_pipeline'],
            'remixt': docker_containers['remixt'],
        },
    }

    return {'count_haps': params}


def get_breakpoint_params(reference_dir, reference, version):
    referencedata = config_reference.get_cluster_reference_data(reference_dir, reference)

    docker_containers = config_reference.containers(version)['docker']

    params = {
        'memory': {'low': 4, 'med': 6, 'high': 16},
        'ref_data_directory': referencedata['destruct_ref_data'],
        'destruct_config': {
            'genome_fasta': referencedata['ref_genome'],
            'genome_fai': referencedata['ref_genome'] + '.fai',
            'gtf_filename': referencedata['destruct_gtf_file'],
        },
        'docker': {
            'single_cell_pipeline': docker_containers['single_cell_pipeline'],
            'destruct': docker_containers['destruct'],
            'lumpy': docker_containers['lumpy'],
            'samtools': docker_containers['samtools'],
        },
    }

    return {'breakpoint_calling': params}


def get_sv_genotyping_params(reference_dir, reference, version):
    docker_containers = config_reference.containers(version)['docker']

    referencedata = config_reference.get_cluster_reference_data(reference_dir, reference)

    params = {
        'ref_genome': referencedata['ref_genome'],
        'memory': {'low': 4, 'med': 6, 'high': 16},
        'docker': {
            'single_cell_pipeline': docker_containers['single_cell_pipeline'],
            'svtyper': docker_containers['svtyper']
        },
    }
    return {'sv_genotyping': params}


def get_singlecell_pipeline_config(config_params, override=None):
    reference = config_params["reference"]
    reference_dir = config_params['refdir']
    version = config_params['version']
    cluster = config_params["cluster"]
    if not reference_dir:
        reference_dir = config_reference.get_reference_dir(cluster)

    params = {}

    params.update(
        get_hmmcopy_params(
            reference_dir, reference, config_params["bin_size"],
            config_params["smoothing_function"], version
        )
    )

    params.update(
        get_align_params(
            reference_dir, reference,
            config_params['aligner'],
            version
        )
    )

    params.update(get_annotation_params(reference_dir, reference, version))

    params.update(get_aneufinder_params(reference_dir, reference, version))

    params.update(get_merge_bams_params(reference_dir, reference, cluster, version))

    params.update(get_split_bam_params(reference_dir, reference, version))

    params.update(get_germline_calling_params(reference_dir, reference, version))

    params.update(get_variant_calling_params(reference_dir, reference, version))

    params.update(get_copy_number_calling_params(
        reference_dir, reference, config_params['copynumber_bin_size'], version)
    )

    params.update(get_infer_haps_params(reference_dir, reference, version))

    params.update(get_count_haps_params(reference_dir, reference, version))

    params.update(get_breakpoint_params(reference_dir, reference, version))

    params.update(get_sv_genotyping_params(reference_dir, reference, version))

    params = override_config(params, override)

    return params
