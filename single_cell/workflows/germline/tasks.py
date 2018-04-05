from biowrappers.components.io.vcf.tasks import index_bcf
import pypeliner


def run_samtools_variant_calling(
        bam_file,
        bai_file,
        ref_genome_fasta_file,
        out_file,
        max_depth=int(1e7),
        min_bqual=0,
        min_depth=0,
        min_mqual=0,
        region=None):

    mpileup_cmd = [
        'samtools',
        'mpileup',
        '-ugf', ref_genome_fasta_file,
        '-Q', min_bqual,
        '-q', min_mqual,
        bam_file
    ]

    if region is not None:
        mpileup_cmd.extend(['-r', region])

    bcf_cmd = [
        'bcftools',
        'call',
        '-vmO', 'z',
        '-o', out_file,
    ]

    cmd = []

    cmd.extend(mpileup_cmd)
    cmd.append('|')
    cmd.extend(bcf_cmd)

    pypeliner.commandline.execute(*cmd)

    index_bcf(out_file)
