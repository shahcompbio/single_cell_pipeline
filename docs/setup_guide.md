###  Download reference data for test datasets
```
wget https://singlecelltestsets.s3.amazonaws.com/refdata.tar.gz
tar -xvf refdata.tar.gz
```


### Align

we recommend starting from a blank slate with a fresh conda install or a new conda environment. 

#### setup conda environment
```
wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh
bash Miniconda3-py37_4.8.3-Linux-x86_64.sh
```
```
conda update -n base -c defaults conda -y   
conda install -c shahcompbio -c bioconda -c conda-forge  single_cell_pipeline_align
```

#### download test dataset:
```
wget https://singlecelltestsets.s3.amazonaws.com/alignment.tar.gz
tar -xvf alignment.tar.gz
```

#### generate inputs.yaml file:
```
SA1090-A96213A-R20-C28:
  column: 28
  condition: B
  fastqs:
    HHCJ7CCXY_5.HGTJJCCXY_8.HYG5LCCXY_6.HYG5LCCXY_7.HYG5LCCXY_5:
      fastq_1: testdata/SA1090-A96213A-R20-C28_1.fastq.gz
      fastq_2: testdata/SA1090-A96213A-R20-C28_2.fastq.gz
      sequencing_center: TEST
      sequencing_instrument: ILLUMINA
  img_col: 45
  index_i5: i5-20
  index_i7: i7-28
  pick_met: C1
  primer_i5: GTATAG
  primer_i7: CTATCT
  row: 20
SA1090-A96213A-R20-C62:
  column: 62
  condition: B
  fastqs:
    HHCJ7CCXY_5.HGTJJCCXY_8.HYG5LCCXY_6.HYG5LCCXY_7.HYG5LCCXY_5:
      fastq_1: testdata/SA1090-A96213A-R20-C62_1.fastq.gz
      fastq_2: testdata/SA1090-A96213A-R20-C62_2.fastq.gz
      sequencing_center: TEST
      sequencing_instrument: ILLUMINA
  img_col: 11
  index_i5: i5-20
  index_i7: i7-62
  pick_met: C1
  primer_i5: GTATAG
  primer_i7: AAGCTA
  row: 20
SA1090-A96213A-R22-C43:
  column: 43
  condition: B
  fastqs:
    HHCJ7CCXY_5.HGTJJCCXY_8.HYG5LCCXY_6.HYG5LCCXY_7.HYG5LCCXY_5:
      fastq_1: testdata/SA1090-A96213A-R22-C43_1.fastq.gz
      fastq_2: testdata/SA1090-A96213A-R22-C43_2.fastq.gz
      sequencing_center: TEST
      sequencing_instrument: ILLUMINA
  img_col: 30
  index_i5: i5-22
  index_i7: i7-43
  pick_met: C2
  primer_i5: GCTGTA
  primer_i7: ATTCCG
  row: 22
```

the testdata path must change to point it to the correct output data directory.

#### launch the pipeline:

```
single_cell alignment --input_yaml inputs.yaml \
--library_id A1234A --maxjobs 4 --nocleanup --sentinel_only \
--submit local --loglevel DEBUG \
--tmpdir temp --pipelinedir pipeline \
--out_dir output --bams_dir bams \
--config_override '{"refdir": "refdata"}'
```

### Hmmcopy

we recommend starting from a blank slate with a fresh conda install or a new conda environment. 

#### setup conda environment
```
wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh
bash Miniconda3-py37_4.8.3-Linux-x86_64.sh
```
```
conda update -n base -c defaults conda -y   
conda install -c conda-forge -c bioconda -c shahcompbio single_cell_pipeline_hmmcopy 
```

#### download test dataset:
```
wget https://singlecelltestsets.s3.amazonaws.com/hmmcopy.tar.gz
tar -xvf hmmcopy.tar.gz
```

#### generate inputs.yaml file:
```
SA1090-A96213A-R20-C28:
  bam: testdata/SA1090-A96213A-R20-C28.bam
  column: 28
  condition: B
  img_col: 45
  index_i5: i5-20
  index_i7: i7-28
  pick_met: C1
  primer_i5: GTATAG
  primer_i7: CTATCT
  row: 20
SA1090-A96213A-R20-C62:
  bam: testdata/SA1090-A96213A-R20-C62.bam
  column: 62
  condition: B
  img_col: 11
  index_i5: i5-20
  index_i7: i7-62
  pick_met: C1
  primer_i5: GTATAG
  primer_i7: AAGCTA
  row: 20
SA1090-A96213A-R22-C43:
  bam: testdata/SA1090-A96213A-R22-C43.bam
  column: 43
  condition: B
  img_col: 30
  index_i5: i5-22
  index_i7: i7-43
  pick_met: C2
  primer_i5: GCTGTA
  primer_i7: ATTCCG
  row: 22
  ```
  the testdata path must change to point it to the correct output data directory.

#### launch the pipeline:

```
 single_cell hmmcopy \
 --input_yaml inputs.yaml \
 --library_id A1234A --maxjobs 4 --nocleanup \
 --sentinel_only --submit local --loglevel DEBUG \
 --tmpdir temp --pipelinedir pipeline --out_dir output \
 --config_override '{"refdir": "refdata", "hmmcopy": {"chromosomes": ["6", "8", "17"]}}'
 ```


### Annotation

we recommend starting from a blank slate with a fresh conda install or a new conda environment. 

#### setup conda environment
```
wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh
bash Miniconda3-py37_4.8.3-Linux-x86_64.sh
```
```
conda update -n base -c defaults conda -y   
conda install -c conda-forge -c bioconda -c shahcompbio single_cell_pipeline_annotation 
```

#### download test dataset:
```
wget https://singlecelltestsets.s3.amazonaws.com/annotation.tar.gz
tar -xvf annotation.tar.gz
```

#### generate inputs.yaml file:
```
hmmcopy_metrics: testdata/A96213A_hmmcopy_metrics.csv.gz
hmmcopy_reads: testdata/A96213A_reads.csv.gz
alignment_metrics: testdata/A96213A_alignment_metrics.csv.gz
gc_metrics: testdata/A96213A_gc_metrics.csv.gz
segs_pdf_tar: testdata/A96213A_segs.tar.gz
```

#### launch the pipeline:

```
 single_cell annotation \
 --input_yaml inputs.yaml \
 --library_id A1234A --maxjobs 4 --nocleanup \
 --sentinel_only --submit local --loglevel DEBUG \
 --tmpdir temp --pipelinedir pipeline --out_dir output \
 --config_override '{"refdir": "refdata", "annotation": {"chromosomes": ["6", "8", "17"]}}'
 ```



### Switching to production runs:

#### Reference data
Before you switch over to production and start running the real datasets, please download the full reference dataset and replace the test dataset from step 1.

```
wget https://singlecelltestsets.s3.amazonaws.com/refdata_full_genome.tar.gz
tar -xvf refdata_full_genome.tar.gz
```

#### Config override
update the config overrides to run the pipeline over the full genome

for Hmmcopy and Annotation, the config override in the launch section should be:
 ```
--config_override '{"refdir": "refdata"}'
```

#### Run with HPC batch submit systems

##### nativespec
We need to figure out the nativespec first. This is a string that specifies the format for job submission.

for instance, the following LSF (for juno cluster at MSKCC) job request
```
bsub -n 1 -W 4:00 -R "rusage[mem=5]span[ptile=1]select[type==CentOS7]"
```

will ask for 1 core, 5 gigs and a runtime of 4 hours

the corresponding pipeline nativespec is
```
-n {ncpus} -W {walltime} -R "rusage[mem={mem}]span[ptile={ncpus}]select[type==CentOS7]"
```


##### launch arguments
please add the following arguments to the launch command

###### Juno cluster (LSF):
```
--submit lsf --nativespec ' -n {ncpus} -W {walltime} -R "rusage[mem={mem}]span[ptile={ncpus}]select[type==CentOS7]"'
```

the pipeline supports *SGE* with `--submit asyncqsub` and *LSF* with `--submit lsf`
