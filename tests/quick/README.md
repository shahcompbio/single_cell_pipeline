# Quick Tests

## Setup

### Downloading test datasets

The most convenient way of accessing the test data is to used `azcopy` version 10 or higher to pull data from azure.  The data should be downloaded to a test_data directory in the single cell pipeline repo to allow compatibility with the inputs.yaml files provided for the test data.

```
cd single_cell_pipeline
azcopy sync https://singlecelltestdata.blob.core.windows.net/data/quick/test_data/ ./test_data/
```

### Create a context_config.yaml

The context config provides information to pypeliner on where to find docker images.  The following should be sufficient:

```
docker:
    server:
    username:
    password:
    mounts:
      refdata: /refdata
      datadrive: /datadrive
      mnt: /mnt
context:
```

## Running the workflows

### Running with docker

Prefix each workflow command with the following to run 

```
docker run -w `pwd` -v `pwd`:`pwd` \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v `which docker`:`which docker` \
    --env PYTHONPATH=`pwd`:`pwd`/pypeliner:`pwd`/biowrappers \
    singlecellpipeline/single_cell_pipeline:v0.5.6 
```

### Alignment workflow

```
single_cell alignment \
    --input_yaml tests/quick/align/inputs.yaml \
    --library_id A96213A \
    --out_dir test_data/align/results \
    --bams_dir test_data/align/bams \
    --submit local --loglevel DEBUG \
    --config_file tests/quick/align/config.yaml \
    --context_config context_config.yaml \
    --tmpdir aligntmp --maxjobs 4 --nocleanup
```

### HMMCopy workflow

```
single_cell hmmcopy \
    --input_yaml tests/quick/hmmcopy/inputs.yaml \
    --library_id A96213A \
    --out_dir test_data/hmmcopy/results \
    --submit local --loglevel DEBUG \
    --config_file tests/quick/hmmcopy/config.yaml \
    --context_config context_config.yaml \
    --tmpdir hmmcopytmp --maxjobs 4 --nocleanup
```

### Annotation workflow

```
single_cell annotation \
    --input_yaml tests/quick/annotation/inputs.yaml \
    --library_id A96213A \
    --out_dir test_data/annotation/results \
    --submit local --loglevel DEBUG \
    --config_file tests/quick/annotation/config.yaml \
    --context_config context_config.yaml \
    --tmpdir annotationtmp --maxjobs 4 --no_corrupt_tree --nocleanup
```

### Split WGS bams workflow

```
single_cell split_wgs_bam \
    --input_yaml tests/quick/split_wgs_bam/inputs.yaml \
    --out_dir test_data/split_wgs_bam/bams \
    --submit local --loglevel DEBUG \
    --config_file tests/quick/split_wgs_bam/config.yaml \
    --context_config context_config.yaml \
    --tmpdir splitwgsbamtmp --maxjobs 4 --nocleanup
```

### Merge cell bams workflow

```
single_cell merge_cell_bams \
    --input_yaml tests/quick/merge_cell_bams/inputs.yaml \
    --out_dir test_data/merge_cell_bams/bams \
    --submit local --loglevel DEBUG \
    --config_file tests/quick/merge_cell_bams/config.yaml \
    --context_config context_config.yaml \
    --tmpdir mergecellbamstmp --maxjobs 4 --nocleanup
```

### Breakpoint calling workflow

```
single_cell breakpoint_calling \
    --input_yaml tests/quick/breakpoint_calling/inputs.yaml \
    --out_dir test_data/breakpoint_calling/results \
    --submit local --loglevel DEBUG \
    --config_file tests/quick/breakpoint_calling/config.yaml \
    --context_config context_config.yaml \
    --tmpdir breakpointcallingtmp --maxjobs 4 --nocleanup
```

### Variant calling workflow

```
single_cell variant_calling \
    --input_yaml tests/quick/variant_calling/inputs.yaml
    --out_dir test_data/variant_calling/results \
    --submit local --loglevel DEBUG \
    --config_file tests/quick/variant_calling/config.yaml \
    --context_config context_config.yaml \
    --tmpdir variantcallingtmp --maxjobs 4 --nocleanup
```
