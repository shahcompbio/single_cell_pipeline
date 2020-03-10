#!/bin/bash
set -e
set -o pipefail

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`
TAG="${TAG}.beta"


mkdir -p SPLIT_WGS_BAM/ref_test_data

docker run -v $PWD:$PWD -w $PWD $3/azurecli:v0.0.1 \
  az storage blob download-batch -s split-bam -d SPLIT_WGS_BAM/ref_test_data --account-name $1 --account-key $2

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v /usr/bin/docker:/usr/bin/docker --rm \
  $3/single_cell_pipeline:$TAG \
  single_cell split_wgs_bam --input_yaml tests/jenkins/split_wgs_bam/inputs.yaml \
  --maxjobs 4 --nocleanup --sentinel_only  \
  --context_config tests/jenkins/context_config.yaml \
  --submit local --loglevel DEBUG \
  --tmpdir SPLIT_WGS_BAM/temp \
  --pipelinedir SPLIT_WGS_BAM/pipeline --submit local --out_dir SPLIT_WGS_BAM/output

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v /usr/bin/docker:/usr/bin/docker --rm \
  $3/single_cell_pipeline:$TAG \
  python tests/jenkins/split_wgs_bam/test_split_wgs_bam.py SPLIT_WGS_BAM/output SPLIT_WGS_BAM/ref_test_data/refdata
