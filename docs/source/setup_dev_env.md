
#### Requirements:

this guide assumes that you have
1. write access to wgspipeline dockerhub org
2. docker installed on the node
3. root access


# Setup


### create working dir and clone code

```
mkdir /devsrc
cd /devsrc
git clone https://github.com/shahcompbio/wgs.git
```

### build docker container

save the following in `dockerfile` file in `/devsrc`
```
# build on top of out base image
FROM wgspipeline/python_wgs:v0.0.1

# Install any needed packages specified in requirements.txt
RUN rm -rf /opt/conda/lib/python2.7/site-packages/pypeliner* /opt/conda/lib/python2.7/site-packages/wgs* /opt/conda/lib/python2.7/site-packages/biowrappers*
RUN pip install git+https://github.com/shahcompbio/pypeliner.git@master
RUN pip install git+https://bitbucket.org/aroth85/biowrappers.git@singlecell
RUN pip install dill
RUN pip install fpdf
RUN pip install intervaltree

#if you plan to shell into container 
RUN apt-get update && apt-get install vim # or editor of your choosing

# Make port 80 available to the world outside this container
EXPOSE 80

# Define environment variable
ENV NAME wgs

ENV PYTHONPATH /devsrc/wgs

# Run app.py when the container launches
CMD ["wgs"]
```

build:
`docker build -t wgspipeline/wgs:v0.0.2 .`

push:
`docker push wgspipeline/wgs:v0.0.2`


### Context config

save the following yaml to a file.


```
docker:
    server: 'docker.io'
    username: null
    password: null
    mounts:
      refdata: /refdata
      datadrive: /devsrc
```

NOTE: you can set the server to `null` if all containers are available locally.

To pull all containers locally, use the `pull_wgs_containers.py` script. The containers are specified in `containers.yaml`.

##### reference data
You need to download the reference data to `/refdata' directory on the node. The data can be downloaded from blob. 

the data is in `wgscomputedev` storage account in `refdata` container.

##### Test Data
Test data for all parts of the wgs pipeline with known ground truths is available in azure storage in the storage account wgstestdata. 
```
wgstestdata\
    alignment-testdata\ #testdata for wgs alignment
        R1.fastq.gz # large test fastqs
        R2.fastq.gz #
        small.R1.fastq.gz # small test fastqs
        small.R2.fastq.gz #
    breakpoint-calling-testdata\ #testdata for wgs breakpoint_calling
        sv_testset\
            bams\ 
                large.bam # large test bam with svs
                medium.bam # medium test bam with svs
                small.bam # small test bam with svs
                normal.bam # no sv test bam
            groundtruth\
                large_bam_breakpoints # expected calls for large bam
                medium_bam_breakpoints.csv # expected calls for medium bam
                small_bam_breakpoints.csv # expected calls for small bam
    cna-testdata\  # testdata for wgs copynumber_calling
        cna\ # cna test data
            config.yaml # config file used for running test bams, pass to wgs copynumber_calling with --config_file argument
            HCC1395BL_chr15_snps.bam # test bam representing normal sample
            HCC1395_chr15_snps.bam # test bam representing tumour sample
            inputs.yaml # input file used for running test bams, pass to wgs copynumber_calling with #input_yaml argument
            targets.tsv # expected output 
        refdata\ # reference data needed to run the cna test set 
            ... # reference files specified in config.yaml
    refdata\
           ... # reference data
    variant-calling-testdata\ #testdata for wgs variant_calling
        data \
            normal.bam
            variants.bam
        configs \
```

#### Alignment pipeline:

1. test data:

The input yaml should have the following format.
```
SA123:
  fastqs:
    L001:
      fastq1: /devsrc/testdata/R1.fastq.gz
      fastq2: /devsrc/testdata/R2.fastq.gz
  bam: results/bams/SA123.bam
```


`R1.fastq.gz` and `R2.fastq.gz` can be found in blob under wgstestdata storage account and alignmenttestdata container. Please contact <grewald@mskcc.org> for access.

2. running

save the following in a shell script:

```
cd /devsrc/wgs && python setup.py develop && cd /devsrc

wgs alignment \
--input_yaml input.yaml \
--out_dir output \
--tmpdir temp \
--pipelinedir pipeline \
--loglevel DEBUG \
--submit local \
--nativespec ' -n {ncpus} -W {walltime} -R "rusage[mem={mem}]span[ptile={ncpus}]select[type==CentOS7]"' \
--sentinel_only  --maxjobs 100 \
--config_override '{"cluster":"juno"}' \
--context_config path_to_context_config_yaml_here
```

and then run the script:

```

docker run -w $PWD -v $PWD:$PWD --rm -v /var/run/docker.sock:/var/run/docker.sock \
-v /usr/bin/docker:/usr/bin/docker  wgspipeline/wgs:v0.0.1 bash run.sh
```

#### Variant Calling Pipeline

1. test data

The input yaml should have the following format
```
SA123:
  normal: /devsrc/testdata/normal.bam
  normal_id: normal
  tumour: /devsrc/testdata/variants.bam
  tumour_id: variants

```

`normal.bam` And `variants.bam` can be found in the `wgstestdata` storage container in the `variant-calling-testdata blob`. The `config.yaml` file used to run the test data can also be found in the `variant-calling-testdata` blob. Please contact <grewald@mskcc.org> for access.

`config.yaml`:

```
variant_calling:
  annotation_params:
    cosmic_params:
      db: /refdata/databases/CosmicMutantExport.sorted.vcf.gz
    dbsnp_params:
      db: /refdata/databases/dbsnp_142.human_9606.all.vcf.gz
    mutation_assessor_params:
      db: /refdata/databases/MA.hg19_v2/
    snpeff_params:
      snpeff_config: /refdata/wgs_pipeline/snpEff.config
    thousandgen_params:
      db: /refdata/databases/1000G_release_20130502_genotypes.vcf.gz
    mappability_ref: /refdata/wgs_pipeline/mask_regions_blacklist_crg_align36_table.txt
  chromosomes:
    - '22'
  docker:
    museqportrait: wgspipeline/museqportrait:v0.0.1
    mutationseq: wgspipeline/mutationseq:v0.0.1
    strelka: wgspipeline/strelka:v0.0.1
    vcftools: wgspipeline/vcftools:v0.0.1
    vizutils: wgspipeline/vizutils:v0.0.1
    wgs: wgspipeline/wgs:v0.0.2
  museq_params:
    baseq_threshold: 20
    buffer_size: 2G
    coverage: 4
    indl_threshold: 0.05
    mapq_threshold: 10
    normal_variant: 25
    purity: 70
    threshold: 0.5
    tumour_variant: 2
    verbose: true
  parse_museq:
    chromosomes:
    - '22'
    keep_1000gen: true
    keep_cosmic: true
    keep_dbsnp: true
    mappability_ref: /refdata/wgs_pipeline/mask_regions_blacklist_crg_align36_table.txt
    pr_threshold: 0.85
    remove_duplicates: false
  parse_strelka:
    chromosomes:
    - '22'
    keep_1000gen: true
    keep_dbsnp: true
    mappability_ref: /refdata/wgs_pipeline/mask_regions_blacklist_crg_align36_table.txt
    remove_duplicates: false
  plot_params:
    dbsnp_params:
      db: /refdata/databases/dbsnp_142.human_9606.all.vcf.gz
    refdata_single_sample: /refdata/wgs_pipeline/single_sample_plot_reference.h5
    thousandgen_params:
      db: /refdata/databases/1000G_release_20130502_genotypes.vcf.gz
    threshold: 0.5
  reference: /refdata/gr37.fasta
  split_size: 10000000.0
```
Note: To use the test data, you must run the pipeline on just chrom `22` and use the file `gr37.fasta` as the reference. 

2. running

save the following in a shell script:

```
cd /devsrc/wgs && python setup.py develop && cd /devsrc

wgs variant_calling \ 
--input_yaml input.yaml 
--out_dir output \
--tmpdir temp \
--pipelinedir pipeline \
--loglevel DEBUG \
--submit local \
--nativespec ' -n {ncpus} -W {walltime} -R "rusage[mem={mem}]span[ptile={ncpus}]select[type==CentOS7]"' \
--sentinel_only --maxjobs 100 \
--config_override '{"cluster":"juno"}' \ 
--config_file config.yaml \
--context_config context.yaml 

```

and then run the script:

`docker run -w $PWD -v $PWD:$PWD --rm -v /var/run/docker.sock:/var/run/docker.sock \
-v /usr/bin/docker:/usr/bin/docker  wgspipeline/wgs:v0.0.1 bash run.sh`

#### Breakpoint Calling

1. test data

The input yaml should have the following format:

```
SA123:
  normal: /devsrc/testdata/normal.bam
  normal_id: normal
  tumour: /devsrc/testdata/small.bam
  tumour_id: breakpoints
```
`normal.bam` And `small.bam` can be found in the `wgstestdata` storage container in the `breakpoint-calling-testdata blob`. The `config.yaml` file used to run the test data can also be found in the `breakpoint-calling-testdata` blob. Please contact <grewald@mskcc.org> for access.

`config.yaml`:

```
sv_calling:
  docker:
    destruct: wgspipeline/destruct:v0.0.4
    lumpy: wgspipeline/lumpy:v0.0.1
    samtools: wgspipeline/samtools:v0.0.2
    vizutils: wgspipeline/vizutils:v0.0.1
    wgs: wgspipeline/wgs:v0.0.2
  extractSplitReads_BwaMem: lumpy_extractSplitReads_BwaMem
  lumpyexpress: lumpyexpress
  mappability_ref: /refdata/wgs_pipeline/mask_regions_blacklist_crg_align36_table.txt
  parse_destruct:
    breakdistance_threshold: 30
    case_id: null
    chromosomes:
    - '22'
    deletion_size_threshold: 1000
    foldback_threshold: 30000
    gene_locations: null
    genes: null
    mappability_ref: /refdata/wgs_pipeline/mask_regions_blacklist_crg_align36_table_destruct.txt
    normal_id: null
    project: null
    readsupport_threshold: 4
    tumour_id: null
    types: null
  parse_lumpy:
    chromosomes:
    - '22'
    confidence_interval_size: 500
    deletion_size_threshold: 0
    foldback_threshold: null
    mappability_ref: null
    normal_id: null
    project: null
    tumour_id: null
    tumour_read_support_threshold: 0
  refdata_destruct: /refdata/reference-destruct/
  samtools: samtoolsb
```

Note: To use the test data, you must run the pipeline using chromosome `22`.

2. running

save the following to a shell script:

```
cd /devsrc/wgs && python setup.py develop && cd /devsrc

wgs breakpoint_calling \ 
--input_yaml input.yaml 
--out_dir output \
--tmpdir temp \
--pipelinedir pipeline \
--loglevel DEBUG \
--submit local \
--nativespec ' -n {ncpus} -W {walltime} -R "rusage[mem={mem}]span[ptile={ncpus}]select[type==CentOS7]"' \
--sentinel_only --maxjobs 100 \
--config_override '{"cluster":"juno"}' \ 
--config_file config.yaml \
--context_config context.yaml 

```


and then run the script:

`docker run -w $PWD -v $PWD:$PWD --rm -v /var/run/docker.sock:/var/run/docker.sock \
-v /usr/bin/docker:/usr/bin/docker  wgspipeline/wgs:v0.0.1 bash run.sh`

#### Copy Number Calling

1. test data

The input yaml should have the following format:

```
HCC1395:
  normal: /devsrc/testdata/HCC1395BL_chr15_snps.bam
  normal_id: HCC1395BL
  tumour: /devsrc/testdata/HCC1395_chr15_snps.bam
  tumour_id: HCC1395
  target_list: /devsrc/testdata/targets.tsv

```
`HCC1395BL_chr15_snps.bam` And `HCC1395_chr15_snps.bam` can be found in the `wgstestdata` storage container in the `cna-testdata blob`. The `config.yaml` file used to run the test data can also be found in the `cna-testdata` blob. Please contact <grewald@mskcc.org> for access.

`config.yaml`:

```
cna_calling:
  chromosomes:
  - '15'
  correction:
    gc: /refdata/wgs/ref_data/GRCh37-lite.gc.ws_1000.wig
  dbsnp_positions: /refdata/databases/common_all_dbSNP138.pos
  docker:
    mutationseq: wgspipeline/mutationseq:v0.0.1
    remixt: wgspipeline/remixt:v0.0.2
    titan: titan:test
    vcftools: wgspipeline/vcftools:v0.0.1
    vizutils: wgspipeline/vizutils:v0.0.1
    wgs: wgs:test
    remixt: amcpherson/remixt:v0.5.7
  min_num_reads: 5
  museq_params:
    baseq_threshold: 10
    buffer_size: 2G
    coverage: 4
    indl_threshold: 0.05
    mapq_threshold: 10
    normal_variant: 25
    purity: 70
    threshold: 0.85
    tumour_variant: 2
    verbose: true
  parse_titan:
    chromosomes:
    - '15'
    genes: null
    segment_size_threshold: 5000
    types: null
  pygenes_gtf: /refdata/ref_data/GRCh37-lite.gtf
  readcounter:
    q: 0
    w: 1000
  reference_genome: /refdata/ref_data/Homo_sapiens.GRCh37.70.dna.chromosomes.fa
  remixt_refdata: /refdata/reference-remixt
  split_size: 10000000.0
  threads: 8
  titan_intervals:
  - num_clusters: 1
    ploidy: 2
  - num_clusters: 2
    ploidy: 2
  titan_params:
    alpha_high: 20000
    alpha_k: 15000
    chrom: 'NULL'
    estimate_ploidy: 'TRUE'
    genome_type: NCBI
    map: /refdata/ref_data/GRCh37-lite.map.ws_1000.wig
    max_copynumber: 8
    max_depth: 1000
    max_iters: 50
    myskew: 0
    normal_estimate_method: map
    normal_param_nzero: 0.5
    num_cores: 4
    pseudo_counts: 1.0e-300
    symmetric: 'TRUE'
    txn_exp_len: 1.0e+16
    txn_z_strength: 1000000.0
    y_threshold: 20
containers:
  docker: &id001
    bwa: wgs/bwa:v0.0.1
    destruct: wgs/destruct:v0.0.2
    lumpy: wgs/lumpy:v0.0.1
    museqportrait: wgs/museqportrait:v0.0.1
    mutationseq: wgs/mutationseq:v0.0.1
    picard: wgs/picard:v0.0.1
    remixt: wgs/remixt:v0.0.2
    samtools: wgs/samtools:v0.0.1
    snpeff: wgs/vcftools:v0.0.1
    strelka: wgs/strelka:v0.0.1
    titan: wgs/titan:v0.0.1
    vcftools: wgs/vcftools:v0.0.1
    vizutils: wgs/vizutils:v0.0.1
    wgs: wgs/wgs:v0.0.2
  singularity: {}
docker_containers:
  docker:
globals:
  memory:
    high: 15
    low: 5
    med: 10
  threads: 8
docker:
```
Note: To use the test data, you must run the pipeline using chromosome `15`.


2. running

save the following to a shell script:

```
cd /devsrc/wgs && python setup.py develop && cd /devsrc

wgs copynumber_calling \ 
--input_yaml input.yaml 
--out_dir output \
--tmpdir temp \
--pipelinedir pipeline \
--loglevel DEBUG \
--submit local \
--nativespec ' -n {ncpus} -W {walltime} -R "rusage[mem={mem}]span[ptile={ncpus}]select[type==CentOS7]"' \
--sentinel_only --maxjobs 100 \
--config_override '{"cluster":"juno"}' \ 
--config_file config.yaml \
--context_config context.yaml 

```


and then run the script:

`docker run -w $PWD -v $PWD:$PWD --rm -v /var/run/docker.sock:/var/run/docker.sock \
-v /usr/bin/docker:/usr/bin/docker  wgspipeline/wgs:v0.0.1 bash run.sh`
