import single_cell


def get_version():
    version = single_cell.__version__
    # strip setuptools metadata
    version = version.split("+")[0]
    return version


def containers():
    version = get_version()
    docker_images = {
        'bwa': 'scp/bwa:v0.0.1', 'samtools': 'scp/samtools:v0.0.2',
        'python_base': 'scp/python_base:v0.0.1', 'picard': 'scp/picard:v0.0.1',
        'single_cell_pipeline': 'scp/single_cell_pipeline:v{}'.format(version),
        'gatk': 'scp/gatk:v0.0.1', 'fastqc': 'scp/fastqc:v0.0.1',
        'hmmcopy': 'scp/hmmcopy:v0.0.2', 'aneufinder': 'scp/aneufinder:v0.0.1',
        'strelka': 'scp/strelka:v0.0.1', 'mutationseq': 'scp/mutationseq:v0.0.2',
        'vcftools': 'scp/vcftools:v0.0.1', 'snpeff': 'scp/vcftools:v0.0.1',
        'titan': 'scp/titan:v0.0.1', 'remixt': 'scp/remixt:v{}'.format(version),
        'destruct': 'scp/destruct:v{}'.format(version), 'trimgalore': 'scp/trimgalore:v0.0.1',
        'lumpy': 'scp/lumpy:v0.0.2',
        'biobloom': 'scp/biobloom:v0.0.2',
    }

    singularity = {}

    return {'docker': docker_images, 'singularity': singularity}


def reference_data_shahlab(reference):
    if reference == 'grch37':
        classifier_training_data = '/shahlab/pipelines/reference/classifier_training_data.h5'
        gc_wig_file = {
            500000: '/shahlab/pipelines/reference/GRCh37-lite.gc.ws_500000.wig',
            1000: '/shahlab/pipelines/reference/GRCh37-lite.gc.ws_1000.wig',
        }
        map_wig_file = {
            500000: '/shahlab/pipelines/reference/GRCh37-lite.map.ws_125_to_500000.wig',
            1000: '/shahlab/pipelines/reference/GRCh37-lite.map.ws_1000.wig',
        }
        exclude_list = '/shahlab/pipelines/reference/repeats.satellite.regions'
        ref_genome = '/shahlab/pipelines/reference/GRCh37-lite.fa'
        ref_genome_index = '/shahlab/pipelines/reference/GRCh37-lite.fa.fai'
        ref_type = reference
        gc_windows = '/shahlab/pipelines/reference/gc_windows.txt'
        copynumber_ref_data = '/shahlab/pipelines/remixt_ref_data_dir'
        chrom_info_filename = '/shahlab/pipelines/remixt_ref_data_dir/chromInfo.txt.gz'
        chromosomes = get_chromosomes('grch37')
        copynumber_ref_data = '/shahlab/pipelines/remixt_ref_data_dir'
        reference_gc_qc = '/shahlab/pipelines/reference/reference_gc_grch37.csv'
        databases = {
            'mappability': {
                'local_path': '/shahlab/pipelines/reference/wgEncodeCrgMapabilityAlign50mer.bigWig',
            },
            'cosmic': {
                'local_path': '/shahlab/pipelines/reference/cosmic_v75.vcf.gz',
            },
            'dbsnp': {
                'local_path': '/shahlab/pipelines/reference/dbsnp_b146_GRCh37p13.vcf.gz',
            }
        }
        biobloom_filters = [
            "/shahlab/pipelines/reference/GCF_002021735.1_Okis_V1_genomic.bf",
            "/shahlab/pipelines/reference/GRCh37-lite.bf",
            "/shahlab/pipelines/reference/mm10_build38_mouse.bf"
        ]
    else:
        classifier_training_data = '/shahlab/pipelines/reference/classifier_training_data.h5'
        gc_wig_file = {
            500000: '//shahlab/pipelines/reference/mm10_build38_mouse.gc.ws_500000.wig',
            1000: '/data/not/available',
        }
        map_wig_file = {
            500000: '//shahlab/pipelines/reference/mm10_build38_mouse.map.rl_125_ws_500000.wig',
            1000: '/data/not/available',
        }
        exclude_list = None
        ref_genome = '/shahlab/pipelines/reference/mm10_build38_mouse.fasta'
        ref_genome_index = '/shahlab/pipelines/reference/mm10_build38_mouse.fasta.fai'
        ref_type = reference
        gc_windows = None
        copynumber_ref_data = '/shahlab/pipelines/remixt_ref_data_dir'
        chrom_info_filename = '/data/not/available'
        chromosomes = get_chromosomes('mm10')
        copynumber_ref_data = '/shahlab/pipelines/remixt_ref_data_dir'
        reference_gc_qc = None
        databases = {
            'mappability': {
                'local_path': None,
            },
            'cosmic': {
                'local_path': None,
            },
            'dbsnp': {
                'local_path': None,
            }
        }
        biobloom_filters = [
             "/shahlab/pipelines/reference/GCF_002021735.1_Okis_V1_genomic.bf",
             "/shahlab/pipelines/reference/GRCh37-lite.bf",
             "/shahlab/pipelines/reference/mm10_build38_mouse.bf"
        ]

    return locals()


def reference_data_azure(reference):
    if reference == "grch37":
        classifier_training_data = '/refdata/classifier_training_data.h5'
        gc_wig_file = {
            500000: '/refdata/GRCh37-lite.gc.ws_500000.wig',
            1000: '/refdata/GRCh37-lite.gc.ws_1000.wig',
        }
        map_wig_file = {
            500000: '/refdata/GRCh37-lite.map.ws_125_to_500000.wig',
            1000: '/refdata/GRCh37-lite.map.ws_1000.wig',
        }
        exclude_list = '/refdata/repeats.satellite.regions'
        ref_genome = '/refdata/GRCh37-lite.fa'
        ref_type = reference
        gc_windows = '/refdata/gc_windows.txt'
        copynumber_ref_data = '/refdata/'
        chrom_info_filename = '/refdata/chromInfo.txt.gz'
        chromosomes = get_chromosomes('grch37')
        copynumber_ref_data = '/refdata/'
        reference_gc_qc = '/refdata/reference_gc_grch37.csv'
        databases = {
            'mappability': {
                'local_path': '/refdata/wgEncodeCrgMapabilityAlign50mer.bigWig',
            },
            'cosmic': {
                'local_path': '/refdata/cosmic_v75.vcf.gz',
            },
            'dbsnp': {
                'local_path': '/refdata/dbsnp_b146_GRCh37p13.vcf.gz',
            }
        }
        biobloom_filters = [
            "/refdata/GCF_002021735.1_Okis_V1_genomic.bf",
            "/refdata/GRCh37-lite.bf",
            "/refdata/mm10_build38_mouse.bf"
        ]

    else:
        classifier_training_data = '/refdata/classifier_training_data.h5'
        gc_wig_file = {
            500000: '/refdata/mm10_build38_mouse.gc.ws_500000.wig',
            1000: '/data/not/available',
        }
        map_wig_file = {
            500000: '/refdata/mm10_build38_mouse.map.rl_125_ws_500000.wig',
            1000: '/data/not/available',
        }
        exclude_list = None
        ref_genome = '/refdata/mm10_build38_mouse.fasta'
        ref_type = reference
        gc_windows = None
        copynumber_ref_data = '/refdata/'
        chrom_info_filename = '/data/not/available'
        chromosomes = get_chromosomes('mm10')
        copynumber_ref_data = '/refdata/remixt_ref_data_dir'
        reference_gc_qc = None
        databases = {
            'mappability': {
                'local_path': None,
            },
            'cosmic': {
                'local_path': None,
            },
            'dbsnp': {
                'local_path': None,
            }
        }
        biobloom_filters = [
            "/refdata/GCF_002021735.1_Okis_V1_genomic.bf",
            "/refdata/GRCh37-lite.bf",
            "/refdata/mm10_build38_mouse.bf"
        ]

    return locals()


def get_chromosomes(reference):
    if reference == 'grch37':
        return map(str, range(1, 23)) + ['X', 'Y']
    elif reference == 'mm10':
        return map(str, range(1, 20)) + ['X', 'Y']
    else:
        raise Exception("unknown reference genome type {}".format(reference))
