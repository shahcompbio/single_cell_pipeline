import logging
import os

logger = logging.getLogger()
logger.setLevel(logging.INFO)
LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

logging.info("Total number of input files: " + str(len(snakemake.input)))


for i in range(len(snakemake.input)):

    logging.info("Sample: " + snakemake.input[i])

    logging.info("Copying vcf file")
    os.system("""cp {inputvcf} results/indelstemp/{pseudobulk_group}/{sample_id}/maf/{sample_id}_{indelid}.vcf.gz
              """.format(indelid = i,
                         inputvcf = snakemake.input[i],
                         pseudobulk_group = snakemake.wildcards.pseudobulk_group,
                         sample_id = snakemake.wildcards.sample_id))

    logging.info("Unzipping vcf file")
    os.system("""gunzip results/indelstemp/{pseudobulk_group}/{sample_id}/maf/{sample_id}_{indelid}.vcf.gz
              """.format(indelid = i,
                         pseudobulk_group = snakemake.wildcards.pseudobulk_group,
                         sample_id = snakemake.wildcards.sample_id))

    logging.info("Running vcf2maf")
    os.system("""vcf2maf.pl --input-vcf results/indelstemp/{pseudobulk_group}/{sample_id}/maf/{sample_id}_{indelid}.vcf \
                 --output-maf results/indelstemp/{pseudobulk_group}/{sample_id}/maf/{sample_id}_{indelid}.maf \
                 --ref-fasta {genomeref} \
                 --filter-vcf 0 --vep-path /opt/vep/src/ensembl-vep  \
                --vep-data {vepdata}""".format(indelid = i,
                                                inputvcf = snakemake.input[i],
                                                pseudobulk_group = snakemake.wildcards.pseudobulk_group,
                                                sample_id = snakemake.wildcards.sample_id,
                                                vepdata = snakemake.params.vepdata,
                                                genomeref = snakemake.params.genomeref))

    logging.info("Remove vcf file")
    os.system("""rm results/indelstemp/{pseudobulk_group}/{sample_id}/maf/{sample_id}_{indelid}.vcf
              """.format(indelid = i,
                         pseudobulk_group = snakemake.wildcards.pseudobulk_group,
                         sample_id = snakemake.wildcards.sample_id))
