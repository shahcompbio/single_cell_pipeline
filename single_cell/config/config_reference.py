import os

from single_cell.utils import helpers


def containers():
    version = helpers.get_version()

    docker_images = {
        'bwa': 'singlecellpipeline/bwa:v0.0.2',
        'samtools': 'singlecellpipeline/samtools:v0.0.3',
        'python_base': 'singlecellpipeline/python_base:v0.1.0',
        'picard': 'singlecellpipeline/picard:v0.0.3',
        'single_cell_pipeline': 'singlecellpipeline/single_cell_pipeline:v{}'.format(version),
        'gatk': 'singlecellpipeline/gatk:v0.0.1',
        'fastqc': 'singlecellpipeline/fastqc:v0.0.2',
        'hmmcopy': 'singlecellpipeline/hmmcopy:v0.0.5',
        'aneufinder': 'singlecellpipeline/aneufinder:v0.0.1',
        'strelka': 'singlecellpipeline/strelka:v0.0.3',
        'mutationseq': 'singlecellpipeline/mutationseq:v0.0.3',
        'vcftools': 'singlecellpipeline/vcftools:v0.0.2',
        'snpeff': 'singlecellpipeline/vcftools:v0.0.2',
        'titan': 'singlecellpipeline/titan:v0.0.1',
        'remixt': 'singlecellpipeline/remixt:v0.0.2',
        'destruct': 'singlecellpipeline/destruct:v0.0.2',
        'trimgalore': 'singlecellpipeline/trimgalore:v0.0.2',
        'lumpy': 'singlecellpipeline/lumpy:v0.0.3',
        'cell_cycle_classifier': 'singlecellpipeline/cell_cycle_classifier:v0.0.1',
        'biobloom': 'singlecellpipeline/biobloom:v0.0.2',
        'corrupt_tree': 'singlecellpipeline/corrupt_tree:v0.0.1',
        'fastq_screen': 'singlecellpipeline/fastq_screen:v0.0.2',
        'svtyper': 'singlecellpipeline/svtyper:v0.0.1'
    }

    return {'docker': docker_images}


def get_reference_data(reference, rootdir):
    # salmon
    if reference == 'GCF_002021735':
        ref_genome = os.path.join(rootdir, 'salmon/GCF_002021735.1_Okis_V1_genomic.fna')
    elif reference == "grch37":
        classifier_training_data = os.path.join(rootdir, 'human/classifier_training_data.h5')
        gc_wig_file = {
            500000: os.path.join(rootdir, 'human/GRCh37-lite.gc.ws_500000.wig'),
            1000: os.path.join(rootdir, 'human/GRCh37-lite.gc.ws_1000.wig'),
        }
        map_wig_file = {
            500000: os.path.join(rootdir, 'human/GRCh37-lite.map.ws_125_to_500000.wig'),
            1000: os.path.join(rootdir, 'human/GRCh37-lite.map.ws_1000.wig'),
        }
        exclude_list = os.path.join(rootdir, 'human/repeats.satellite.regions')
        ref_genome = os.path.join(rootdir, 'human/GRCh37-lite.fa')
        gc_windows = os.path.join(rootdir, 'human/gc_windows.txt')
        copynumber_ref_data = os.path.join(rootdir, 'human/')
        chrom_info_filename = os.path.join(rootdir, 'human/chromInfo.txt.gz')
        chromosomes = get_chromosomes('grch37')
        destruct_ref_data = os.path.join(rootdir, 'human/')
        destruct_gtf_file = os.path.join(rootdir, 'human/GRCh37-lite.gtf')
        reference_gc_qc = os.path.join(rootdir, 'human/reference_gc_grch37.csv')
        databases = {
            'mappability': {
                'local_path': os.path.join(rootdir, 'human/wgEncodeCrgMapabilityAlign50mer.bigWig'),
            },
            'cosmic': {
                'local_path': os.path.join(rootdir, 'human/cosmic_v75.vcf.gz'),
            },
            'dbsnp': {
                'local_path': os.path.join(rootdir, 'human/dbsnp_b146_GRCh37p13.vcf.gz'),
            }
        }

    else:
        classifier_training_data = os.path.join(rootdir, 'mouse/classifier_training_data.h5')
        gc_wig_file = {
            500000: os.path.join(rootdir, 'mouse/mm10_build38_mouse.gc.ws_500000.wig'),
            1000: '/data/not/available',
        }
        map_wig_file = {
            500000: os.path.join(rootdir, 'mouse/mm10_build38_mouse.map.rl_125_ws_500000.wig'),
            1000: '/data/not/available',
        }
        exclude_list = None
        ref_genome = os.path.join(rootdir, 'mouse/mm10_build38_mouse.fasta')
        gc_windows = None
        copynumber_ref_data = os.path.join(rootdir, 'mouse/')
        chrom_info_filename = '/data/not/available'
        chromosomes = get_chromosomes('mm10')
        destruct_ref_data = os.path.join(rootdir, 'mouse/')
        destruct_gtf_file = os.path.join(rootdir, 'mouse/mm10_build38_mouse.gtf')
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

    return locals()


def get_reference_dir(cluster):
    if cluster == "shahlab":
        reference_dir = '/shahlab/pipelines/reference/singlecellpipeline'
    elif cluster == "juno":
        reference_dir = '/juno/work/shah/reference/singlecellpipeline'
    else:
        reference_dir = '/refdata'

    return reference_dir


def get_cluster_reference_data(refdir, reference):
    return get_reference_data(reference, refdir)


def get_chromosomes(reference):
    if reference == 'grch37':
        return [str(val) for val in range(1, 23)] + ['X', 'Y']
    elif reference == 'mm10':
        return [str(val) for val in range(1, 20)] + ['X', 'Y']
    else:
        raise Exception("unknown reference genome type {}".format(reference))
