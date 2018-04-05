from pypeliner.workflow import Workflow

import csv
import pypeliner
import pysam
import vcf_tasks
import tasks
import os


def create_snv_allele_counts_for_vcf_targets_workflow(
        bam_file,
        vcf_file,
        out_file,
        chromosomes=default_chromosomes,
        count_duplicates=False,
        hdf5_output=True,
        min_bqual=0,
        min_mqual=0,
        split_size=int(1e7),
        table_name='snv_allele_counts',
        vcf_to_bam_chrom_map=None):

    if hdf5_output:
        merged_file = mgd.File(out_file)

    else:
        merged_file = mgd.TempFile('merged.h5')

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='get_regions',
        ret=mgd.TempOutputObj('regions_obj', 'regions'),
        func=utils.get_vcf_regions,
        args=(
            mgd.InputFile(vcf_file),
            split_size,
        ),
        kwargs={
            'chromosomes': chromosomes,
        },
    )

    workflow.transform(
        name='get_snv_allele_counts_for_vcf_targets',
        axes=('regions',),
        ctx=med_ctx,
        func=tasks.get_snv_allele_counts_for_vcf_targets,
        args=(
            mgd.InputFile(bam_file),
            mgd.InputFile(vcf_file),
            mgd.TempOutputFile('counts.h5', 'regions'),
            table_name
        ),
        kwargs={
            'count_duplicates': count_duplicates,
            'min_bqual': min_bqual,
            'min_mqual': min_mqual,
            'region': mgd.TempInputObj('regions_obj', 'regions'),
            'vcf_to_bam_chrom_map': vcf_to_bam_chrom_map,
        }
    )

    workflow.transform(
        name='merge_snv_allele_counts',
        ctx=med_ctx,
        func=hdf5_tasks.concatenate_tables,
        args=(
            mgd.TempInputFile('counts.h5', 'regions'),
            merged_file.as_output(),
        ),
        kwargs={
            'in_memory': False,
        }
    )

    if not hdf5_output:
        workflow.transform(
            name='convert_to_tsv',
            ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2},
            func=hdf5_tasks.convert_hdf5_to_tsv,
            args=(
                merged_file.as_input(),
                table_name,
                mgd.OutputFile(out_file),
            ),
            kwargs={
                'compress': True,
            }
        )

    return workflow


def create_germline_workflow(
        normal_bam_files,
        normal_bai_files,
        tumour_bam_files,
        tumour_bai_files,
        ref_genome_fasta_file,
        vcf_file,
        parsed_csv,
        config,
        regions,
        chromosomes=default_chromosomes,
        split_size=int(1e7),
        use_depth_thresholds=True):

    workflow = Workflow()

    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('chrom'),
        value=chromosomes,
    )
    
    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('regions'),
        value=regions,
    )
              
    workflow.transform(
        name='count_fasta_bases',
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 },
        func=tasks.count_fasta_bases,
        args=(
            ref_genome_fasta_file,
            pypeliner.managed.TempOutputFile('ref_base_counts.tsv')
        )
    )

    workflow.transform(
        name="get_chrom_sizes",
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 },
        func=tasks.get_known_chromosome_sizes,
        ret=pypeliner.managed.TempOutputObj('known_sizes'),
        args=(
              pypeliner.managed.TempInputFile('ref_base_counts.tsv'),
              chromosomes
        )
    )

    workflow.transform(
        name='run_samtools_variant_calling',
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2, 'ncpus':1, 'pool_id': config['pools']['standard']},
        axes=('regions',),
        func=tasks.run_samtools_variant_calling,
        args=(
            pypeliner.managed.InputFile('normal.bam', 'regions', fnames=normal_bam_files),
            pypeliner.managed.InputFile('normal.bam.bai', 'regions', fnames=normal_bai_files),
            pypeliner.managed.InputFile(ref_genome_fasta_file),
            pypeliner.managed.TempOutputFile('variants.vcf.gz', 'regions'),
        ),
        kwargs={
            'region': pypeliner.managed.InputInstance('regions'),
        },
    )

    workflow.transform(
        name='concatenate_variants',
        ctx={'mem': 2},
        func=vcf_tasks.concatenate_vcf,
        args=(
            pypeliner.managed.TempInputFile('variants.vcf.gz', 'regions'),
            pypeliner.managed.OutputFile(vcf_file),
        ),
    )

    workflow.subworkflow(
        name='read_counts',
        axes=('sample_id',),
        func=biowrappers.components.variant_calling.snv_allele_counts.create_snv_allele_counts_for_vcf_targets_workflow,
        args=(
            mgd.InputFile('bam', 'sample_id', fnames=bam_filenames),
            mgd.InputFile(germline_vcf_filename),
            mgd.OutputFile(counts_template, 'sample_id'),
        ),
        kwargs={
            'table_name': mgd.Template('/counts/{sample_id}', 'sample_id'),
        },
    )
    









    workflow.transform(
        name='add_indel_filters',
        axes=('chrom',),
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 },
        func=tasks.filter_indel_file_list,
        args=(
            pypeliner.managed.TempInputFile('somatic.indels.unfiltered.vcf', 'regions'),
            pypeliner.managed.TempInputFile('strelka.stats', 'regions'),
            pypeliner.managed.TempInputFile('somatic.indels.unfiltered.vcf.window', 'regions'),
            pypeliner.managed.TempOutputFile('somatic.indels.filtered.vcf', 'chrom'),
            pypeliner.managed.InputInstance("chrom"),
            pypeliner.managed.TempInputObj('known_sizes'),
            regions
        ),
        kwargs={'use_depth_filter': use_depth_thresholds}
    )
  
    workflow.transform(
        name='add_snv_filters',
        axes=('chrom',),
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 },
        func=tasks.filter_snv_file_list,
        args=(
            pypeliner.managed.TempInputFile('somatic.snvs.unfiltered.vcf', 'interval', axes_origin=[]),
            pypeliner.managed.TempInputFile('strelka.stats', 'interval', axes_origin=[]),
            pypeliner.managed.TempOutputFile('somatic.snvs.filtered.vcf', 'chrom'),
            pypeliner.managed.InputInstance("chrom"),
            pypeliner.managed.TempInputObj('known_sizes'),
            regions,
        ),
        kwargs={'use_depth_filter': use_depth_thresholds}
    )
    
    workflow.transform(
        name='merge_indels',
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 },
        func=vcf_tasks.concatenate_vcf,
        args=(
            pypeliner.managed.TempInputFile('somatic.indels.filtered.vcf', 'chrom'),
            pypeliner.managed.TempOutputFile('somatic.indels.filtered.vcf.gz')
        )
    )
    
    workflow.transform(
        name='merge_snvs',
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 },
        func=vcf_tasks.concatenate_vcf,
        args=(
            pypeliner.managed.TempInputFile('somatic.snvs.filtered.vcf', 'chrom'),
            pypeliner.managed.TempOutputFile('somatic.snvs.filtered.vcf.gz')
        )
    )
    
    workflow.transform(
        name='filter_indels',
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 },
        func=vcf_tasks.filter_vcf,
        args=(
            pypeliner.managed.TempInputFile('somatic.indels.filtered.vcf.gz'),
            pypeliner.managed.TempOutputFile('somatic.indels.passed.vcf')
        )
    )
    
    workflow.transform(
        name='filter_snvs',
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 },
        func=vcf_tasks.filter_vcf,
        args=(
            pypeliner.managed.TempInputFile('somatic.snvs.filtered.vcf.gz'),
            pypeliner.managed.TempOutputFile('somatic.snvs.passed.vcf')
        )
    )
    
    workflow.transform(
        name='finalise_indels',
        ctx={'pool_id': config['pools']['standard'], 'ncpus':1},
        func=vcf_tasks.finalise_vcf,
        args=(
            pypeliner.managed.TempInputFile('somatic.indels.passed.vcf'),
            pypeliner.managed.OutputFile(indel_vcf_file)
        )
    )
    
    workflow.transform(
        name='finalise_snvs',
        ctx={'pool_id': config['pools']['standard'], 'ncpus':1},
        func=vcf_tasks.finalise_vcf,
        args=(
            pypeliner.managed.TempInputFile('somatic.snvs.passed.vcf'),
            pypeliner.managed.OutputFile(snv_vcf_file)
        )
    )
    
    workflow.transform(
        name='parse_strelka_snv',
        ctx={'mem': 10, 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.parse_strelka,
        args=(
              pypeliner.managed.InputFile(snv_vcf_file),
              pypeliner.managed.OutputFile(parsed_snv_csv),
              )
        )
    
    workflow.transform(
        name='parse_strelka_indel',
        ctx={'mem': 10, 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.parse_strelka,
        args=(
              pypeliner.managed.InputFile(indel_vcf_file),
              pypeliner.managed.OutputFile(parsed_indel_csv),
              )
         )

    return workflow


def get_chromosomes(bam_file, chromosomes=None):

    chromosomes = _get_chromosomes(bam_file, chromosomes)

    return dict(zip(chromosomes, chromosomes))


def _get_chromosomes(bam_file, chromosomes=None):
    bam = pysam.Samfile(bam_file, 'rb')

    if chromosomes is None:
        chromosomes = bam.references

    else:
        chromosomes = chromosomes

    return [str(x) for x in chromosomes]


def get_coords(bam_file, chrom, split_size):

    coords = {}

    bam = pysam.Samfile(bam_file, 'rb')

    chrom_lengths = dict(zip(bam.references, bam.lengths))

    length = chrom_lengths[chrom]

    lside_interval = range(1, length + 1, split_size)

    rside_interval = range(split_size, length + split_size, split_size)

    for coord_index, (beg, end) in enumerate(zip(lside_interval, rside_interval)):
        coords[coord_index] = (beg, end)

    return coords


def get_known_chromosome_sizes(bam_file, bai_file, size_file, chromosomes):
    chromosomes = _get_chromosomes(bam_file, chromosomes)

    sizes = {}

    with open(size_file, 'r') as fh:
        reader = csv.DictReader(fh, ['path', 'chrom', 'known_size', 'size'], delimiter='\t')

        for row in reader:
            if row['chrom'] not in chromosomes:
                continue

            sizes[row['chrom']] = int(row['known_size'])

    return sizes
