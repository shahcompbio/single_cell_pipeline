'''
Created on Jun 6, 2018

@author: dgrewal
'''

def get_chromosomes(reference):
    if reference == 'grch37':
        return map(str, range(1,23)) + ['X','Y']
    elif reference == 'mm10':
        return map(str, range(1,20)) + ['X','Y']
    else:
        raise Exception("unknown reference genome type {}".format(reference))


def get_reference_azure():

    common_options = {
        'max_cores': None,
        'mutationseq': '/usr/local/museq/',
        'mutationseq_model': '/usr/local/museq/models_anaconda/model_v4.1.2_anaconda_sk_0.13.1.npz',
        'mutationseq_python': '/usr/local/miniconda2/envs/museq/bin/python',
        'one_split_job': True,

    }

    human_data = {
        'grch37': {
            'gc_wig_file': {
                500000: '/refdata/GRCh37-lite.gc.ws_500000.wig',
                1000: '/refdata/GRCh37-lite.gc.ws_1000.wig',
            },
            'map_wig_file': {
                500000: '/refdata/GRCh37-lite.map.ws_125_to_500000.wig',
                1000: '/refdata/GRCh37-lite.map.ws_1000.wig',
            },
            'exclude_list': '/refdata/repeats.satellite.regions',
            'ref_genome': '/refdata/GRCh37-lite.fa',
            'ref_genome_index': '/refdata/GRCh37-lite.fa.fai',
            'gc_windows': '/refdata/gc_windows.txt',
            'copynumber_ref_data': '/refdata/',
            'chrom_info_filename': '/refdata/chromInfo.txt.gz',
            'chromosomes': get_chromosomes('grch37'),
        },
    }

    mouse_data = {
        'mm10': {
            'gc_wig_file': {
                500000: '/refdata/mm10_build38_mouse.gc.ws_500000.wig',
                1000: '/data/not/available',
            },
            'map_wig_file': {
                500000: '/refdata/mm10_build38_mouse.map.rl_125_ws_500000.wig',
                1000: '/data/not/available',
            },
            'exclude_list': None,
            'ref_genome': '/refdata/mm10_build38_mouse.fasta',
            'ref_genome_index': '/refdata/mm10_build38_mouse.fasta.fai',
            'gc_windows': None,
            'copynumber_ref_data': '/refdata/',
            'chrom_info_filename': '/data/not/available',
            'chromosomes': get_chromosomes('mm10'),
        },
    }

    data = {}
    data.update(common_options)
    data.update(human_data)
    data.update(mouse_data)

    return {'azure': data}


def get_reference_shahlab():

    common_options = {
        'max_cores': 8,
        'mutationseq': '/shahlab/pipelines/apps_centos6/mutationSeq_4.3.7_anaconda/',
        'mutationseq_model': '/shahlab/pipelines/apps_centos6/mutationSeq_4.3.7_anaconda//model_v4.1.2_anaconda_sk_0.13.1.npz',
        'mutationseq_python': 'python',
        'one_split_job': False,
    }

    human_data = {
        'grch37': {
            'gc_wig_file': {
                500000: '/shahlab/pipelines/reference/GRCh37-lite.gc.ws_500000.wig',
                1000: '/shahlab/pipelines/reference/GRCh37-lite.gc.ws_1000.wig',
            },
            'map_wig_file': {
                500000: '/shahlab/pipelines/reference/GRCh37-lite.map.ws_125_to_500000.wig',
                1000: '/shahlab/pipelines/reference/GRCh37-lite.map.ws_1000.wig',
            },
            'exclude_list': '/shahlab/pipelines/reference/repeats.satellite.regions',
            'ref_genome': '/shahlab/pipelines/reference/GRCh37-lite.fa',
            'ref_genome_index': '/shahlab/pipelines/reference/GRCh37-lite.fa.fai',
            'gc_windows': '/shahlab/pipelines/reference/gc_windows.txt',
            'copynumber_ref_data': '/shahlab/pipelines/remixt_ref_data_dir',
            'chrom_info_filename': '/shahlab/pipelines/remixt_ref_data_dir/chromInfo.txt.gz',
            'chromosomes': get_chromosomes('grch37'),
        },
    }

    mouse_data = {
        'mm10': {
            'gc_wig_file': {
                500000: '//shahlab/pipelines/reference/mm10_build38_mouse.gc.ws_500000.wig',
                1000: '/data/not/available',
            },
            'map_wig_file': {
                500000: '//shahlab/pipelines/reference/mm10_build38_mouse.map.rl_125_ws_500000.wig',
                1000: '/data/not/available',
            },
            'exclude_list': None,
            'ref_genome': '//shahlab/pipelines/reference/mm10_build38_mouse.fasta',
            'ref_genome_index': '//shahlab/pipelines/reference/mm10_build38_mouse.fasta.fai',
            'gc_windows': None,
            'copynumber_ref_data': '/shahlab/pipelines/remixt_ref_data_dir',
            'chrom_info_filename': '/data/not/available',
            'chromosomes': get_chromosomes('mm10'),
        },
    }

    data = {}
    data.update(common_options)
    data.update(human_data)
    data.update(mouse_data)

    return {'shahlab': data}


def get_databases():

    cosmic_db = {
        'download_method': 'sftp',
        'user_name': 'awm3@sfu.ca',
        'password': 'shahlabith',
        'host': 'sftp-cancer.sanger.ac.uk',
        'remote_paths': {
            'coding': '/files/grch37/cosmic/v75/VCF/CosmicCodingMuts.vcf.gz',
            'non_coding': '/files/grch37/cosmic/v75/VCF/CosmicNonCodingVariants.vcf.gz',
        },
        'local_path': '/refdata/cosmic_v75.vcf.gz',
    }

    dbsnp_db = {
        'url': 'ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b146_GRCh37p13/VCF/common_all_20151104.vcf.gz',
        'local_path': '/refdata/dbsnp_b146_GRCh37p13.vcf.gz',
    }

    mapp_db = {
        'url': 'http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/release3/wgEncodeCrgMapabilityAlign50mer.bigWig',
        'local_path': '/refdata/wgEncodeCrgMapabilityAlign50mer.bigWig',
    }

    reference = {
        'url': 'http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa',
        'local_path': '/refdata/GRCh37-lite.fa',
    }

    snpeff = {"db": 'GRCh37.75'}

    databases = {
        'cosmic': cosmic_db,
        'dbsnp': dbsnp_db,
        'mappability': mapp_db,
        'ref_genome': reference,
        'snpeff': snpeff
    }

    databases = {
        "databases": databases,
    }

    return databases

def get_const():

    params = {
        "min_bqual":20,
        "min_mqual":20,
        "count_unpaired": False,        
        "normal_contamination": [0.2, 0.4, 0.6, 0.8],
        'num_clusters': [1, 2],
        'ploidy': [1, 2, 3, 4],
        'multipliers': [1, 2, 3, 4, 5, 6],
        'map_cutoff': 0.9,
        'e': 0.999999,
        'eta': 50000,
        'g': 3,
        'lambda': 20,
        'nu': 2.1,
        'num_states': 12,
        's': 1,
        'strength': 1000,
        'kappa': '100,100,700,100,25,25,25,25,25,25,25,25',
        'm': '0,1,2,3,4,5,6,7,8,9,10,11',
        'mu': '0,1,2,3,4,5,6,7,8,9,10,11',
        'good_cells': [['median_hmmcopy_reads_per_bin', 'ge', 50]],
        'memory': {'low': 2, 'med': 6, 'high': 18},
        'adapter': 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC',
        'adapter2': 'CTGTCTCTTATACACATCTGACGCTGCCGACGA',
        'split_size': 10000000,

    }
    
    return params

def reference_data():
    data = {}
    data.update(get_reference_azure())
    data.update(get_reference_shahlab())
    data.update(get_databases())
    data.update(get_const())
    return data


def extract_from_reference(keys):

    data = reference_data()

    keydata = data
    for key in keys:
        keydata = keydata[key]

    return keydata
