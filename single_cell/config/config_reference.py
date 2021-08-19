import os


def get_reference_data(reference, rootdir):
    # salmon
    if reference == 'GCF_002021735':
        ref_genome = os.path.join(rootdir, 'salmon/GCF_002021735.1_Okis_V1_genomic.fna')
    elif reference == "grch37":
        vep = {
            'reference_dir': os.path.join(rootdir, 'human', 'vep'),
            'reference_fasta': os.path.join(
                rootdir, 'human', 'vep', 'homo_sapiens', '99_GRCh37',
                'Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz'
            ),
            'reference_filter_vcf': os.path.join(
                rootdir, 'human', 'vep', 'ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz'
            )
        }
        classifier_training_data = os.path.join(rootdir, 'human/classifier_training_data.h5')
        fastqscreen_training_data = os.path.join(rootdir, 'human/fastqscreen_training_data.csv')
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
        qc_gtf_file = os.path.join(rootdir, 'human/Homo_sapiens.GRCh37.73.gtf')
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
            },
            'snpeff': {
                'local_path': os.path.join(rootdir, 'snpeff/data/'),
                'db': 'GRCh37.75'
            }
        }

    elif reference == 'mm10':
        vep = {
            'reference_dir': None,
            'reference_fasta': None,
            'reference_filter_vcf': None
        }
        classifier_training_data = os.path.join(rootdir, 'mouse/classifier_training_data.h5')
        fastqscreen_training_data = os.path.join(rootdir, 'human/fastqscreen_training_data.csv')
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
        qc_gtf_file = None
        reference_gc_qc = None
        qc_gtf_file = None
        databases = {
            'mappability': {
                'local_path': os.path.join(rootdir, 'mouse/k50.Umap.MultiTrackMappability.bw'),
            },
            'mgp_indel': {
                'local_path': os.path.join(rootdir, 'mouse/mgp.v3.indels.rsIDdbSNPv137.vcf.gz'),
            },
            'mgp_snv': {
                'local_path': os.path.join(rootdir, 'mouse/mgp.v3.snps.rsIDdbSNPv137.vcf.gz'),
            },
            'snpeff': {
                'local_path': os.path.join(rootdir, 'snpeff/data'),
                'db': 'mm10'
            }
        }
    else:
        raise Exception('unknown reference')

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
