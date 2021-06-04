
### 1. Start a new directory for Single Cell Pipeline

```
mkdir single_cell_pipeline
cd  single_cell_pipeline
```

_note down absolute path to this directory. We'll need this later._



###  2. Download reference data for test datasets
```
wget https://singlecelltestsets.s3.amazonaws.com/refdata.tar.gz
tar -xvf refdata.tar.gz
```
_note down absolute path to the extracted reference directory. We'll need this later._


### 3. version
Track down the pipeline version. Please refer to [CHANGELOG](../../CHANGELOG.md) for release notes and versions.


### 4. Alignment

create a separate directory for alignment:
```
mkdir ALIGN && cd ALIGN
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
      trim: true
  img_col: 45
  index_i5: i5-20
  index_i7: i7-28
  pick_met: C1
  primer_i5: GTATAG
  primer_i7: CTATCT
  row: 20
  sample_id: SA1090
  library_id: A96213A
SA1090-A96213A-R20-C62:
  column: 62
  condition: B
  fastqs:
    HHCJ7CCXY_5.HGTJJCCXY_8.HYG5LCCXY_6.HYG5LCCXY_7.HYG5LCCXY_5:
      fastq_1: testdata/SA1090-A96213A-R20-C62_1.fastq.gz
      fastq_2: testdata/SA1090-A96213A-R20-C62_2.fastq.gz
      sequencing_center: TEST
      trim: true
  img_col: 11
  index_i5: i5-20
  index_i7: i7-62
  pick_met: C1
  primer_i5: GTATAG
  primer_i7: AAGCTA
  row: 20
  sample_id: SA1090
  library_id: A96213A
SA1090-A96213A-R22-C43:
  column: 43
  condition: B
  fastqs:
    HHCJ7CCXY_5.HGTJJCCXY_8.HYG5LCCXY_6.HYG5LCCXY_7.HYG5LCCXY_5:
      fastq_1: testdata/SA1090-A96213A-R22-C43_1.fastq.gz
      fastq_2: testdata/SA1090-A96213A-R22-C43_2.fastq.gz
      sequencing_center: TEST
      trim: true
  img_col: 30
  index_i5: i5-22
  index_i7: i7-43
  pick_met: C2
  primer_i5: GCTGTA
  primer_i7: ATTCCG
  row: 22
  sample_id: SA1090
  library_id: A96213A
```

the testdata path must change to point it to the correct output data directory.

#### launch the pipeline:

create `runner.sh`
```
single_cell alignment --input_yaml inputs.yaml \
--library_id A1234A --maxjobs 4 --nocleanup --sentinel_only \
--submit local --loglevel DEBUG \
--tmpdir temp --pipelinedir pipeline \
--out_dir output --bams_dir bams \
--config_override '{"refdir": "<REFERENCE_DIR>"}' \
--sequencing_center TEST
```

please replace `<REFERENCE_DIR>` with path to the reference data we extracted in step 2.


##### Singularity:

build the singularity container
```
singularity build scp_alignment.sif docker://quay.io/singlecellpipeline/single_cell_pipeline_alignment:v<version>
```
please replace `<version>` with the version identified in step 3.

launch:
```
singularity run --bind <MOUNT_DIR>  scp_alignment.sif sh runner.sh
```
please replace `<MOUNT_DIR>` with path to the directory in step 1.


##### Docker:

```
 docker run -w $PWD -v <MOUNT_DIR>:<MOUNT_DIR> quay.io/singlecellpipeline/single_cell_pipeline_alignment:v<version> sh runner.sh
```
please replace `<version>` with the version identified in step 3.
please replace `<MOUNT_DIR>` with path to the directory in step 1.


### 4. Hmmcopy

#### setup conda environment

create a separate directory for hmmcopy:
```
mkdir HMMCOPY && cd HMMCOPY
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
  sample_id: SA1090
  library_id: A96213A
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
  sample_id: SA1090
  library_id: A96213A
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
  sample_id: SA1090
  library_id: A96213A
  ```
  the testdata path must change to point it to the correct output data directory.

#### launch the pipeline:

create `runner.sh`
```
 single_cell hmmcopy \
 --input_yaml inputs.yaml \
 --library_id A1234A --maxjobs 4 --nocleanup \
 --sentinel_only --submit local --loglevel DEBUG \
 --tmpdir temp --pipelinedir pipeline --out_dir output \
 --config_override '{"refdir": "<REFERENCE_DIR>", "hmmcopy": {"chromosomes": ["6", "8", "17"]}}'
```

please replace `<REFERENCE_DIR>` with path to the reference data we extracted in step 2.


##### Singularity:

build the singularity container
```
singularity build scp_hmmcopy.sif docker://quay.io/singlecellpipeline/single_cell_pipeline_hmmcopy:v<version>
```
please replace `<version>` with the version identified in step 3.


launch:
```
singularity run --bind <MOUNT_DIR>  scp_hmmcopy.sif sh runner.sh
```
please replace `<MOUNT_DIR>` with path to the directory in step 1.


##### Docker:

```
 docker run -w $PWD -v <MOUNT_DIR>:<MOUNT_DIR> quay.io/singlecellpipeline/single_cell_pipeline_hmmcopy:v<version> sh runner.sh
```
please replace `<version>` with the version identified in step 3.
please replace `<MOUNT_DIR>` with path to the directory in step 1.



### 5. Annotation


create a separate directory for annotation:
```
mkdir ANNOTATION && cd ANNOTATION
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

create `runner.sh`
```
 single_cell annotation \
 --input_yaml inputs.yaml \
 --library_id A1234A --maxjobs 4 --nocleanup \
 --sentinel_only --submit local --loglevel DEBUG \
 --tmpdir temp --pipelinedir pipeline --out_dir output \
 --config_override '{"refdir": "<REFERENCE_DIR>", "annotation": {"chromosomes": ["6", "8", "17"]}}'
```

please replace `<REFERENCE_DIR>` with path to the reference data we extracted in step 2.


##### Singularity:

build the singularity container
```
singularity build scp_annotation.sif docker://quay.io/singlecellpipeline/single_cell_pipeline_annotation:v<version>
```
please replace `<version>` with the version identified in step 3.


launch:
```
singularity run --bind <MOUNT_DIR>  scp_annotation.sif sh runner.sh
```
please replace `<MOUNT_DIR>` with path to the directory in step 1.


##### Docker:

```
 docker run -w $PWD -v <MOUNT_DIR>:<MOUNT_DIR> quay.io/singlecellpipeline/single_cell_pipeline_annotation:v<version> sh runner.sh
```
please replace `<version>` with the version identified in step 3.
please replace `<MOUNT_DIR>` with path to the directory in step 1.



### 6. merge cells


create a separate directory for annotation:
```
mkdir MERGE_CELLS && cd MERGE_CELLS
```

#### download test dataset:
```
wget https://singlecelltestsets.s3.amazonaws.com/merge_cells.tar.gz
tar -xvf merge_cells.tar.gz
```

#### generate inputs.yaml file:
```
cell_bams:
  SA1090-A96213A-R20-C28:
    bam: merge_cells/SA1090-A96213A-R20-C28.bam
  SA1090-A96213A-R20-C62:
    bam: merge_cells/SA1090-A96213A-R20-C62.bam
  SA1090-A96213A-R22-C43:
    bam: merge_cells/SA1090-A96213A-R22-C43.bam
  SA1090-A96213A-R22-C44:
    bam: merge_cells/SA1090-A96213A-R22-C44.bam
  SA1090-A96213A-R24-C12:
    bam: merge_cells/SA1090-A96213A-R24-C12.bam
  SA1090-A96213A-R24-C20:
    bam: merge_cells/SA1090-A96213A-R24-C20.bam
  SA1090-A96213A-R24-C58:
    bam: merge_cells/SA1090-A96213A-R24-C58.bam
  SA1090-A96213A-R25-C14:
    bam: merge_cells/SA1090-A96213A-R25-C14.bam
  SA1090-A96213A-R25-C22:
    bam: merge_cells/SA1090-A96213A-R25-C22.bam
  SA1090-A96213A-R25-C40:
    bam: merge_cells/SA1090-A96213A-R25-C40.bam
  SA1090-A96213A-R25-C64:
    bam: merge_cells/SA1090-A96213A-R25-C64.bam
  SA1090-A96213A-R26-C49:
    bam: merge_cells/SA1090-A96213A-R26-C49.bam
  SA1090-A96213A-R26-C50:
    bam: merge_cells/SA1090-A96213A-R26-C50.bam
  SA1090-A96213A-R26-C61:
    bam: merge_cells/SA1090-A96213A-R26-C61.bam
  SA1090-A96213A-R26-C64:
    bam: merge_cells/SA1090-A96213A-R26-C64.bam
  SA1090-A96213A-R27-C14:
    bam: merge_cells/SA1090-A96213A-R27-C14.bam
  SA1090-A96213A-R27-C21:
    bam: merge_cells/SA1090-A96213A-R27-C21.bam
  SA1090-A96213A-R27-C45:
    bam: merge_cells/SA1090-A96213A-R27-C45.bam
  SA1090-A96213A-R28-C23:
    bam: merge_cells/SA1090-A96213A-R28-C23.bam
  SA1090-A96213A-R28-C39:
    bam: merge_cells/SA1090-A96213A-R28-C39.bam
  SA1090-A96213A-R28-C64:
    bam: merge_cells/SA1090-A96213A-R28-C64.bam
  SA1090-A96213A-R29-C18:
    bam: merge_cells/SA1090-A96213A-R29-C18.bam
  SA1090-A96213A-R29-C47:
    bam: merge_cells/SA1090-A96213A-R29-C47.bam
  SA1090-A96213A-R29-C59:
    bam: merge_cells/SA1090-A96213A-R29-C59.bam
  SA1090-A96213A-R29-C62:
    bam: merge_cells/SA1090-A96213A-R29-C62.bam
  SA1090-A96213A-R30-C15:
    bam: merge_cells/SA1090-A96213A-R30-C15.bam
  SA1090-A96213A-R30-C35:
    bam: merge_cells/SA1090-A96213A-R30-C35.bam
  SA1090-A96213A-R30-C55:
    bam: merge_cells/SA1090-A96213A-R30-C55.bam
  SA1090-A96213A-R31-C09:
    bam: merge_cells/SA1090-A96213A-R31-C09.bam
  SA1090-A96213A-R31-C37:
    bam: merge_cells/SA1090-A96213A-R31-C37.bam
  SA1090-A96213A-R32-C39:
    bam: merge_cells/SA1090-A96213A-R32-C39.bam
  SA1090-A96213A-R32-C41:
    bam: merge_cells/SA1090-A96213A-R32-C41.bam
  SA1090-A96213A-R32-C65:
    bam: merge_cells/SA1090-A96213A-R32-C65.bam
  SA1090-A96213A-R32-C66:
    bam: merge_cells/SA1090-A96213A-R32-C66.bam
  SA1090-A96213A-R33-C31:
    bam: merge_cells/SA1090-A96213A-R33-C31.bam
  SA1090-A96213A-R33-C38:
    bam: merge_cells/SA1090-A96213A-R33-C38.bam
  SA1090-A96213A-R33-C66:
    bam: merge_cells/SA1090-A96213A-R33-C66.bam
  SA1090-A96213A-R35-C18:
    bam: merge_cells/SA1090-A96213A-R35-C18.bam
  SA1090-A96213A-R35-C24:
    bam: merge_cells/SA1090-A96213A-R35-C24.bam
  SA1090-A96213A-R35-C38:
    bam: merge_cells/SA1090-A96213A-R35-C38.bam
  SA1090-A96213A-R35-C47:
    bam: merge_cells/SA1090-A96213A-R35-C47.bam
```
    

#### launch the pipeline:

create `runner.sh`
```
 single_cell merge_cell_bams \
 --input_yaml inputs.yaml \
 --maxjobs 4 --nocleanup \
 --sentinel_only --submit local --loglevel DEBUG \
 --tmpdir temp --pipelinedir pipeline --out_dir output \
 --config_override '{"refdir": "<REFERENCE_DIR>", "merge_bams":{"chromosomes": ["6", "8", "17"]}}'
```

please replace `<REFERENCE_DIR>` with path to the reference data we extracted in step 2.


##### Singularity:

build the singularity container
```
singularity build scp_merge_cells.sif docker://quay.io/singlecellpipeline/single_cell_pipeline_alignment:v<version>
```
please replace `<version>` with the version identified in step 3.


launch:
```
singularity run --bind <MOUNT_DIR>  scp_merge_cells.sif sh runner.sh
```
please replace `<MOUNT_DIR>` with path to the directory in step 1.


##### Docker:

```
 docker run -w $PWD -v <MOUNT_DIR>:<MOUNT_DIR> quay.io/singlecellpipeline/single_cell_pipeline_alignment:v<version> sh runner.sh
```
please replace `<version>` with the version identified in step 3.
please replace `<MOUNT_DIR>` with path to the directory in step 1.




### 7. split wgs bam


create a separate directory for annotation:
```
mkdir SPLIT_WGS_BAM && cd SPLIT_WGS_BAM
```

#### download test dataset:
```
wget https://singlecelltestsets.s3.amazonaws.com/split_bam.tar.gz
tar -xvf split_bam.tar.gz
```

#### generate inputs.yaml file:
```
normal:
  bam: split_bam/DAH370N_A41086.bam
```
    

#### launch the pipeline:

create `runner.sh`
```
 single_cell split_wgs_bam \
 --input_yaml inputs.yaml \
 --maxjobs 4 --nocleanup \
 --sentinel_only --submit local --loglevel DEBUG \
 --tmpdir temp --pipelinedir pipeline --out_dir output \
 --config_override '{"refdir": "<REFERENCE_DIR>"}'
```

please replace `<REFERENCE_DIR>` with path to the reference data we extracted in step 2.

##### Singularity:

build the singularity container
```
singularity build scp_split_wgs_bam.sif docker://quay.io/singlecellpipeline/single_cell_pipeline_alignment:v<version>
```
please replace `<version>` with the version identified in step 3.


launch:
```
singularity run --bind <MOUNT_DIR>  scp_split_wgs_bam.sif sh runner.sh
```
please replace `<MOUNT_DIR>` with path to the directory in step 1.


##### Docker:

```
 docker run -w $PWD -v <MOUNT_DIR>:<MOUNT_DIR> quay.io/singlecellpipeline/single_cell_pipeline_annotation:v<version> sh runner.sh
```
please replace `<version>` with the version identified in step 3.
please replace `<MOUNT_DIR>` with path to the directory in step 1.


### 8. variant calling


create a separate directory for annotation:
```
mkdir VARIANT_CALLING && cd VARIANT_CALLING
```

#### download test dataset:
```
wget https://singlecelltestsets.s3.amazonaws.com/variant_calling.tar.gz
tar -xvf variant_calling.tar.gz
```

#### generate inputs.yaml file:
```
normal:
  17-1-10000000:
    bam: variant_calling/normal/17-1-10000000.bam
  17-10000001-20000000:
    bam: variant_calling/normal/17-10000001-20000000.bam
  17-20000001-30000000:
    bam: variant_calling/normal/17-20000001-30000000.bam
  17-30000001-40000000:
    bam: variant_calling/normal/17-30000001-40000000.bam
  17-40000001-50000000:
    bam: variant_calling/normal/17-40000001-50000000.bam
  17-50000001-60000000:
    bam: variant_calling/normal/17-50000001-60000000.bam
  17-60000001-70000000:
    bam: variant_calling/normal/17-60000001-70000000.bam
  17-70000001-80000000:
    bam: variant_calling/normal/17-70000001-80000000.bam
  17-80000001-81195210:
    bam: variant_calling/normal/17-80000001-81195210.bam
  6-1-10000000:
    bam: variant_calling/normal/6-1-10000000.bam
  6-100000001-110000000:
    bam: variant_calling/normal/6-100000001-110000000.bam
  6-10000001-20000000:
    bam: variant_calling/normal/6-10000001-20000000.bam
  6-110000001-120000000:
    bam: variant_calling/normal/6-110000001-120000000.bam
  6-120000001-130000000:
    bam: variant_calling/normal/6-120000001-130000000.bam
  6-130000001-140000000:
    bam: variant_calling/normal/6-130000001-140000000.bam
  6-140000001-150000000:
    bam: variant_calling/normal/6-140000001-150000000.bam
  6-150000001-160000000:
    bam: variant_calling/normal/6-150000001-160000000.bam
  6-160000001-170000000:
    bam: variant_calling/normal/6-160000001-170000000.bam
  6-170000001-171115067:
    bam: variant_calling/normal/6-170000001-171115067.bam
  6-20000001-30000000:
    bam: variant_calling/normal/6-20000001-30000000.bam
  6-30000001-40000000:
    bam: variant_calling/normal/6-30000001-40000000.bam
  6-40000001-50000000:
    bam: variant_calling/normal/6-40000001-50000000.bam
  6-50000001-60000000:
    bam: variant_calling/normal/6-50000001-60000000.bam
  6-60000001-70000000:
    bam: variant_calling/normal/6-60000001-70000000.bam
  6-70000001-80000000:
    bam: variant_calling/normal/6-70000001-80000000.bam
  6-80000001-90000000:
    bam: variant_calling/normal/6-80000001-90000000.bam
  6-90000001-100000000:
    bam: variant_calling/normal/6-90000001-100000000.bam
  8-1-10000000:
    bam: variant_calling/normal/8-1-10000000.bam
  8-100000001-110000000:
    bam: variant_calling/normal/8-100000001-110000000.bam
  8-10000001-20000000:
    bam: variant_calling/normal/8-10000001-20000000.bam
  8-110000001-120000000:
    bam: variant_calling/normal/8-110000001-120000000.bam
  8-120000001-130000000:
    bam: variant_calling/normal/8-120000001-130000000.bam
  8-130000001-140000000:
    bam: variant_calling/normal/8-130000001-140000000.bam
  8-140000001-146364022:
    bam: variant_calling/normal/8-140000001-146364022.bam
  8-20000001-30000000:
    bam: variant_calling/normal/8-20000001-30000000.bam
  8-30000001-40000000:
    bam: variant_calling/normal/8-30000001-40000000.bam
  8-40000001-50000000:
    bam: variant_calling/normal/8-40000001-50000000.bam
  8-50000001-60000000:
    bam: variant_calling/normal/8-50000001-60000000.bam
  8-60000001-70000000:
    bam: variant_calling/normal/8-60000001-70000000.bam
  8-70000001-80000000:
    bam: variant_calling/normal/8-70000001-80000000.bam
  8-80000001-90000000:
    bam: variant_calling/normal/8-80000001-90000000.bam
  8-90000001-100000000:
    bam: variant_calling/normal/8-90000001-100000000.bam
tumour:
  17-1-10000000:
    bam: variant_calling/tumour/17-1-10000000.bam
  17-10000001-20000000:
    bam: variant_calling/tumour/17-10000001-20000000.bam
  17-20000001-30000000:
    bam: variant_calling/tumour/17-20000001-30000000.bam
  17-30000001-40000000:
    bam: variant_calling/tumour/17-30000001-40000000.bam
  17-40000001-50000000:
    bam: variant_calling/tumour/17-40000001-50000000.bam
  17-50000001-60000000:
    bam: variant_calling/tumour/17-50000001-60000000.bam
  17-60000001-70000000:
    bam: variant_calling/tumour/17-60000001-70000000.bam
  17-70000001-80000000:
    bam: variant_calling/tumour/17-70000001-80000000.bam
  17-80000001-81195210:
    bam: variant_calling/tumour/17-80000001-81195210.bam
  6-1-10000000:
    bam: variant_calling/tumour/6-1-10000000.bam
  6-100000001-110000000:
    bam: variant_calling/tumour/6-100000001-110000000.bam
  6-10000001-20000000:
    bam: variant_calling/tumour/6-10000001-20000000.bam
  6-110000001-120000000:
    bam: variant_calling/tumour/6-110000001-120000000.bam
  6-120000001-130000000:
    bam: variant_calling/tumour/6-120000001-130000000.bam
  6-130000001-140000000:
    bam: variant_calling/tumour/6-130000001-140000000.bam
  6-140000001-150000000:
    bam: variant_calling/tumour/6-140000001-150000000.bam
  6-150000001-160000000:
    bam: variant_calling/tumour/6-150000001-160000000.bam
  6-160000001-170000000:
    bam: variant_calling/tumour/6-160000001-170000000.bam
  6-170000001-171115067:
    bam: variant_calling/tumour/6-170000001-171115067.bam
  6-20000001-30000000:
    bam: variant_calling/tumour/6-20000001-30000000.bam
  6-30000001-40000000:
    bam: variant_calling/tumour/6-30000001-40000000.bam
  6-40000001-50000000:
    bam: variant_calling/tumour/6-40000001-50000000.bam
  6-50000001-60000000:
    bam: variant_calling/tumour/6-50000001-60000000.bam
  6-60000001-70000000:
    bam: variant_calling/tumour/6-60000001-70000000.bam
  6-70000001-80000000:
    bam: variant_calling/tumour/6-70000001-80000000.bam
  6-80000001-90000000:
    bam: variant_calling/tumour/6-80000001-90000000.bam
  6-90000001-100000000:
    bam: variant_calling/tumour/6-90000001-100000000.bam
  8-1-10000000:
    bam: variant_calling/tumour/8-1-10000000.bam
  8-100000001-110000000:
    bam: variant_calling/tumour/8-100000001-110000000.bam
  8-10000001-20000000:
    bam: variant_calling/tumour/8-10000001-20000000.bam
  8-110000001-120000000:
    bam: variant_calling/tumour/8-110000001-120000000.bam
  8-120000001-130000000:
    bam: variant_calling/tumour/8-120000001-130000000.bam
  8-130000001-140000000:
    bam: variant_calling/tumour/8-130000001-140000000.bam
  8-140000001-146364022:
    bam: variant_calling/tumour/8-140000001-146364022.bam
  8-20000001-30000000:
    bam: variant_calling/tumour/8-20000001-30000000.bam
  8-30000001-40000000:
    bam: variant_calling/tumour/8-30000001-40000000.bam
  8-40000001-50000000:
    bam: variant_calling/tumour/8-40000001-50000000.bam
  8-50000001-60000000:
    bam: variant_calling/tumour/8-50000001-60000000.bam
  8-60000001-70000000:
    bam: variant_calling/tumour/8-60000001-70000000.bam
  8-70000001-80000000:
    bam: variant_calling/tumour/8-70000001-80000000.bam
  8-80000001-90000000:
    bam: variant_calling/tumour/8-80000001-90000000.bam
  8-90000001-100000000:
    bam: variant_calling/tumour/8-90000001-100000000.bam
``` 

#### launch the pipeline:

create `runner.sh`
```
 single_cell variant_calling \
 --input_yaml inputs.yaml \
 --maxjobs 4 --nocleanup \
 --sentinel_only --submit local --loglevel DEBUG \
 --tmpdir temp --pipelinedir pipeline --out_dir output \
 --config_override '{"refdir": "<REFERENCE_DIR>", "variant_calling": {"chromosomes": ["6", "8", "17"]}}'
```

please replace `<REFERENCE_DIR>` with path to the reference data we extracted in step 2.


##### Singularity:

build the singularity container
```
singularity build scp_variant.sif docker://quay.io/singlecellpipeline/single_cell_pipeline_variant:v<version>
```
please replace `<version>` with the version identified in step 3.


launch:
```
singularity run --bind <MOUNT_DIR>  scp_variant.sif sh runner.sh
```
please replace `<MOUNT_DIR>` with path to the directory in step 1.


##### Docker:

```
 docker run -w $PWD -v <MOUNT_DIR>:<MOUNT_DIR> quay.io/singlecellpipeline/single_cell_pipeline_variant:v<version> sh runner.sh
```
please replace `<version>` with the version identified in step 3.
please replace `<MOUNT_DIR>` with path to the directory in step 1.






### 9. breakpoint calling


create a separate directory for annotation:
```
mkdir BREAKPOINT_CALLING && cd BREAKPOINT_CALLING
```

#### download test dataset:
```
wget https://singlecelltestsets.s3.amazonaws.com/breakpoint_calling.tar.gz
tar -xvf breakpoint_calling.tar.gz
```

#### generate inputs.yaml file:
```
normal:
  bam: breakpoint_calling/normal/DAH370N_A41086.bam
tumour:
  SA1090-A96213A-R20-C28:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R20-C28.bam
  SA1090-A96213A-R20-C62:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R20-C62.bam
  SA1090-A96213A-R22-C43:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R22-C43.bam
  SA1090-A96213A-R22-C44:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R22-C44.bam
  SA1090-A96213A-R24-C12:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R24-C12.bam
  SA1090-A96213A-R24-C19:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R24-C19.bam
  SA1090-A96213A-R24-C20:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R24-C20.bam
  SA1090-A96213A-R24-C58:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R24-C58.bam
  SA1090-A96213A-R25-C13:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R25-C13.bam
  SA1090-A96213A-R25-C14:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R25-C14.bam
  SA1090-A96213A-R25-C22:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R25-C22.bam
  SA1090-A96213A-R25-C40:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R25-C40.bam
  SA1090-A96213A-R25-C64:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R25-C64.bam
  SA1090-A96213A-R26-C49:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R26-C49.bam
  SA1090-A96213A-R26-C50:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R26-C50.bam
  SA1090-A96213A-R26-C61:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R26-C61.bam
  SA1090-A96213A-R26-C64:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R26-C64.bam
  SA1090-A96213A-R27-C10:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R27-C10.bam
  SA1090-A96213A-R27-C14:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R27-C14.bam
  SA1090-A96213A-R27-C21:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R27-C21.bam
  SA1090-A96213A-R27-C45:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R27-C45.bam
  SA1090-A96213A-R28-C12:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R28-C12.bam
  SA1090-A96213A-R28-C20:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R28-C20.bam
  SA1090-A96213A-R28-C23:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R28-C23.bam
  SA1090-A96213A-R28-C39:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R28-C39.bam
  SA1090-A96213A-R28-C55:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R28-C55.bam
  SA1090-A96213A-R28-C64:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R28-C64.bam
  SA1090-A96213A-R29-C18:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R29-C18.bam
  SA1090-A96213A-R29-C47:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R29-C47.bam
  SA1090-A96213A-R29-C59:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R29-C59.bam
  SA1090-A96213A-R29-C62:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R29-C62.bam
  SA1090-A96213A-R29-C68:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R29-C68.bam
  SA1090-A96213A-R30-C14:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R30-C14.bam
  SA1090-A96213A-R30-C15:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R30-C15.bam
  SA1090-A96213A-R30-C35:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R30-C35.bam
  SA1090-A96213A-R30-C44:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R30-C44.bam
  SA1090-A96213A-R30-C55:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R30-C55.bam
  SA1090-A96213A-R31-C09:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R31-C09.bam
  SA1090-A96213A-R31-C37:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R31-C37.bam
  SA1090-A96213A-R32-C39:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R32-C39.bam
  SA1090-A96213A-R32-C41:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R32-C41.bam
  SA1090-A96213A-R32-C65:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R32-C65.bam
  SA1090-A96213A-R32-C66:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R32-C66.bam
  SA1090-A96213A-R33-C31:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R33-C31.bam
  SA1090-A96213A-R33-C38:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R33-C38.bam
  SA1090-A96213A-R33-C66:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R33-C66.bam
  SA1090-A96213A-R35-C18:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R35-C18.bam
  SA1090-A96213A-R35-C24:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R35-C24.bam
  SA1090-A96213A-R35-C25:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R35-C25.bam
  SA1090-A96213A-R35-C28:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R35-C28.bam
  SA1090-A96213A-R35-C38:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R35-C38.bam
  SA1090-A96213A-R35-C47:
    bam: breakpoint_calling/tumour/SA1090-A96213A-R35-C47.bam
``` 

#### launch the pipeline:

create `runner.sh`
```
 single_cell breakpoint_calling \
 --input_yaml inputs.yaml \
 --maxjobs 4 --nocleanup \
 --sentinel_only --submit local --loglevel DEBUG \
 --tmpdir temp --pipelinedir pipeline --out_dir output \
 --config_override '{"refdir": "<REFERENCE_DIR>"}'
```

please replace `<REFERENCE_DIR>` with path to the reference data we extracted in step 2.


##### Singularity:

build the singularity container
```
singularity build scp_breakpoint_calling.sif docker://quay.io/singlecellpipeline/single_cell_pipeline_breakpoint:v<version>
```
please replace `<version>` with the version identified in step 3.


launch:
```
singularity run --bind <MOUNT_DIR>  scp_breakpoint_calling.sif sh runner.sh
```
please replace `<MOUNT_DIR>` with path to the directory in step 1.


##### Docker:

```
 docker run -w $PWD -v <MOUNT_DIR>:<MOUNT_DIR> quay.io/singlecellpipeline/single_cell_pipeline_breakpoint:v<version> sh runner.sh
```
please replace `<version>` with the version identified in step 3.
please replace `<MOUNT_DIR>` with path to the directory in step 1.




### 10. infer haplotypes


create a separate directory for annotation:
```
mkdir INFER_HAPS && cd INFER_HAPS
```

#### download test dataset:
```
wget https://singlecelltestsets.s3.amazonaws.com/infer_haps.tar.gz
tar -xvf infer_haps.tar.gz
```

#### generate inputs.yaml file:
```
normal:
    bam: infer_haps/HCC1395BL_chr15.bam
``` 

#### launch the pipeline:

create `runner.sh`
```
 single_cell infer_haps \
 --input_yaml inputs.yaml \
 --maxjobs 4 --nocleanup \
 --sentinel_only --submit local --loglevel DEBUG \
 --tmpdir temp --pipelinedir pipeline --out_dir output \
 --config_override '{"refdir": "<REFERENCE_DIR>", "infer_haps":{"chromosomes":["15"]}}'
```

please replace `<REFERENCE_DIR>` with path to the reference data we extracted in step 2.

##### Singularity:

build the singularity container
```
singularity build scp_infer_haps.sif docker://quay.io/singlecellpipeline/single_cell_pipeline_haplotypes:v<version>
```
please replace `<version>` with the version identified in step 3.


launch:
```
singularity run --bind <MOUNT_DIR>  scp_infer_haps.sif sh runner.sh
```
please replace `<MOUNT_DIR>` with path to the directory in step 1.


##### Docker:

```
 docker run -w $PWD -v <MOUNT_DIR>:<MOUNT_DIR> quay.io/singlecellpipeline/single_cell_pipeline_haplotypes:v<version> sh runner.sh
```
please replace `<version>` with the version identified in step 3.
please replace `<MOUNT_DIR>` with path to the directory in step 1.



### 11. count haplotypes


create a separate directory for annotation:
```
mkdir COUNT_HAPS && cd COUNT_HAPS
```


#### download test dataset:
```
wget https://singlecelltestsets.s3.amazonaws.com/count_haps.tar.gz
tar -xvf count_haps.tar.gz
```

#### generate inputs.yaml file:
```
haplotypes: count_haps/haps.csv.gz
tumour:
  SA607_3X10XB02284-A108843A-R03-C03:
    bam: count_haps/SA607_3X10XB02284-A108843A-R03-C03.bam
  SA607_3X10XB02284-A108843A-R03-C10:
    bam: count_haps/SA607_3X10XB02284-A108843A-R03-C10.bam
  SA607_3X10XB02284-A108843A-R03-C08:
    bam: count_haps/SA607_3X10XB02284-A108843A-R03-C08.bam
  SA607_3X10XB02284-A108843A-R03-C09:
    bam: count_haps/SA607_3X10XB02284-A108843A-R03-C09.bam
``` 

#### launch the pipeline:

create `runner.sh`
```
 single_cell count_haps \
 --input_yaml inputs.yaml \
 --maxjobs 4 --nocleanup \
 --sentinel_only --submit local --loglevel DEBUG \
 --tmpdir temp --pipelinedir pipeline --out_dir output \
 --config_override '{"refdir": "<REFERENCE_DIR>"}'
```

please replace `<REFERENCE_DIR>` with path to the reference data we extracted in step 2.

##### Singularity:

build the singularity container
```
singularity build scp_count_haps.sif docker://quay.io/singlecellpipeline/single_cell_pipeline_haplotypes:v<version>
```
please replace `<version>` with the version identified in step 3.


launch:
```
singularity run --bind <MOUNT_DIR>  scp_count_haps.sif sh runner.sh
```
please replace `<MOUNT_DIR>` with path to the directory in step 1.


##### Docker:

```
 docker run -w $PWD -v <MOUNT_DIR>:<MOUNT_DIR> quay.io/singlecellpipeline/single_cell_pipeline_haplotypes:v<version> sh runner.sh
```
please replace `<version>` with the version identified in step 3.
please replace `<MOUNT_DIR>` with path to the directory in step 1.



### 12. snv genotyping


create a separate directory for annotation:
```
mkdir SNV_GENOTYPING && cd SNV_GENOTYPING
```


#### download test dataset:
```
wget https://singlecelltestsets.s3.amazonaws.com/snv_genotyping.tar.gz
tar -xvf snv_genotyping.tar.gz
```

#### generate inputs.yaml file:
```
vcf_files:
  - snv_genotyping/vcf/museq.vcf.gz
  - snv_genotyping/vcf/strelka_snv.vcf.gz
tumour_cells:
  SA1090:
    A96213A:
      SA1090-A96213A-R20-C28:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R20-C28.bam
      SA1090-A96213A-R22-C43:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R22-C43.bam
      SA1090-A96213A-R22-C44:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R22-C44.bam
      SA1090-A96213A-R24-C12:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R24-C12.bam
      SA1090-A96213A-R24-C20:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R24-C20.bam
      SA1090-A96213A-R24-C58:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R24-C58.bam
      SA1090-A96213A-R25-C14:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R25-C14.bam
      SA1090-A96213A-R25-C22:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R25-C22.bam
      SA1090-A96213A-R25-C40:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R25-C40.bam
      SA1090-A96213A-R25-C64:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R25-C64.bam
      SA1090-A96213A-R26-C49:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R26-C49.bam
      SA1090-A96213A-R26-C50:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R26-C50.bam
      SA1090-A96213A-R26-C64:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R26-C64.bam
      SA1090-A96213A-R27-C14:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R27-C14.bam
      SA1090-A96213A-R27-C21:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R27-C21.bam
      SA1090-A96213A-R27-C45:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R27-C45.bam
      SA1090-A96213A-R28-C23:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R28-C23.bam
      SA1090-A96213A-R28-C39:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R28-C39.bam
      SA1090-A96213A-R28-C64:
        bam: snv_genotyping/cell_bams/SA1090-A96213A-R28-C64.bam
```
#### launch the pipeline:

create `runner.sh`
```
 single_cell snv_genotyping \
 --input_yaml inputs.yaml \
 --maxjobs 4 --nocleanup \
 --sentinel_only --submit local --loglevel DEBUG \
 --tmpdir temp --pipelinedir pipeline --out_dir output \
 --config_override '{"refdir": "<REFERENCE_DIR>"}'
```

please replace `<REFERENCE_DIR>` with path to the reference data we extracted in step 2.

##### Singularity:

build the singularity container
```
singularity build scp_snv_genotyping.sif docker://quay.io/singlecellpipeline/single_cell_pipeline_variant:v<version>
```
please replace `<version>` with the version identified in step 3.


launch:
```
singularity run --bind <MOUNT_DIR>  scp_snv_genotyping.sif sh runner.sh
```
please replace `<MOUNT_DIR>` with path to the directory in step 1.


##### Docker:

```
 docker run -w $PWD -v <MOUNT_DIR>:<MOUNT_DIR> quay.io/singlecellpipeline/single_cell_pipeline_variant:v<version> sh runner.sh
```
please replace `<version>` with the version identified in step 3.
please replace `<MOUNT_DIR>` with path to the directory in step 1.




### 13. Switching to production runs:

#### Reference data
Before you switch over to production and start running the real datasets, please download the full reference dataset and replace the test dataset from step 1.

```
wget https://singlecelltestsets.s3.amazonaws.com/refdata_full_genome.tar.gz
tar -xvf refdata_full_genome.tar.gz
```

update the config overrides to run the pipeline over the full genome, the config override in the launch section should point to the full reference dir:
 ```
--config_override '{"refdir": "refdata"}'
```