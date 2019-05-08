import pypeliner
import pypeliner.managed as mgd
from pypeliner.workflow import Workflow

default_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']


def create_strelka_workflow(
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        indel_vcf_file,
        snv_vcf_file,
        config,
        chromosomes=default_chromosomes
):
    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1,
           'num_retry': 3, 'docker_image': config['docker']['single_cell_pipeline']}

    regions = list(normal_bam_file.keys())
    assert set(tumour_bam_file.keys()) == set(regions)

    regions = [val for val in regions if val.startswith('1-')]

    workflow = Workflow(ctx=ctx)

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
        ctx=dict(mem=2),
        func="single_cell.workflows.strelka.tasks.count_fasta_bases",
        args=(
            ref_genome_fasta_file,
            pypeliner.managed.TempOutputFile('ref_base_counts.tsv'),
        ),
        kwargs={'docker_image': config['docker']['strelka']}
    )

    workflow.transform(
        name="get_chrom_sizes",
        ctx=dict(mem=2),
        func="single_cell.workflows.strelka.tasks.get_known_chromosome_sizes",
        ret=pypeliner.managed.TempOutputObj('known_sizes'),
        args=(
            pypeliner.managed.TempInputFile('ref_base_counts.tsv'),
            chromosomes
        )
    )

    workflow.transform(
        name='get_chromosome_depths',
        axes=('regions',),
        func="single_cell.workflows.strelka.tasks.get_chromosome_depth",
        args=(
            mgd.InputInstance('regions'),
            pypeliner.managed.InputFile("normal.split.bam", "regions", fnames=normal_bam_file, extensions=['.bai']),
            ref_genome_fasta_file,
            mgd.TempOutputFile('chrom_depth.txt', 'regions'),
        ),
        kwargs={'docker_image': config['docker']['strelka']},
    )

    workflow.transform(
        name='merge_chromosome_depths',
        func="single_cell.workflows.strelka.tasks.merge_chromosome_depths",
        args=(
            mgd.TempInputFile('chrom_depth.txt', 'regions'),
            mgd.TempOutputFile('merged_chrom_depth.txt')
        )
    )

    workflow.transform(
        name='call_genome_segment',
        axes=('regions',),
        func="single_cell.workflows.strelka.tasks.call_genome_segment",
        args=(
            mgd.TempInputFile('merged_chrom_depth.txt'),
            pypeliner.managed.InputFile("normal.split.bam", "regions", fnames=normal_bam_file, extensions=['.bai']),
            pypeliner.managed.InputFile("merged_bam", "regions", fnames=tumour_bam_file, extensions=['.bai']),
            ref_genome_fasta_file,
            mgd.TempOutputFile('indels.vcf', 'regions'),
            mgd.TempOutputFile('snvs.vcf', 'regions'),
            mgd.TempSpace('call_genome_segment_tmp', 'regions'),
            mgd.InputInstance('regions'),
            mgd.TempInputObj('known_sizes'),
        ),
        kwargs={
            'is_exome': False,
            'docker_image': config['docker']['strelka']
        }
    )

    workflow.transform(
        name='merge_indels',
        func='single_cell.utils.vcfutils.concatenate_vcf',
        args=(
            mgd.TempInputFile('indels.vcf', 'regions'),
            mgd.TempOutputFile('indels.vcf.gz', extensions=['.tbi', '.csi']),
            mgd.TempSpace("indels_merge")
        ),
        kwargs={'docker_image': config['docker']['vcftools']}
    )

    workflow.transform(
        name='merge_snvs',
        func='single_cell.utils.vcfutils.concatenate_vcf',
        args=(
            mgd.TempInputFile('snvs.vcf', 'regions'),
            mgd.TempOutputFile('snvs.vcf.gz', extensions=['.tbi', '.csi']),
            mgd.TempSpace("snvs_merge")
        ),
        kwargs={'docker_image': config['docker']['vcftools']}
    )

    workflow.transform(
        name='filter_vcf_indel',
        func='single_cell.utils.vcfutils.filter_vcf',
        args=(
            mgd.TempInputFile('indels.vcf.gz', extensions=['.tbi', '.csi']),
            mgd.OutputFile(indel_vcf_file, extensions=['.tbi', '.csi']),
        ),
        kwargs={'docker_image': config['docker']['vcftools']}
    )

    workflow.transform(
        name='filter_vcf_snv',
        func='single_cell.utils.vcfutils.filter_vcf',
        args=(
            mgd.TempInputFile('snvs.vcf.gz', extensions=['.tbi', '.csi']),
            mgd.OutputFile(snv_vcf_file, extensions=['.tbi', '.csi']),
        ),
        kwargs={'docker_image': config['docker']['vcftools']}
    )

    return workflow
