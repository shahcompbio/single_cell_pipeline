#!/bin/bash
set -e
set -o pipefail

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`
DOCKER=`which docker`

mkdir -p SPLIT_WGS_BAM/ref_test_data

docker run -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY -e AWS_DEFAULT_REGION -v $PWD:$PWD -w $PWD $1/awscli:v0.0.1 \
  aws s3 cp s3://singlecelltestsets/TESTDATA_CODEBUILD/split-bam SPLIT_WGS_BAM/ref_test_data --recursive


docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v $DOCKER:$DOCKER --rm \
  $1/single_cell_pipeline:$TAG \
  single_cell split_wgs_bam --input_yaml single_cell/tests/jenkins/split_wgs_bam/inputs.yaml \
  --maxjobs 4 --nocleanup --sentinel_only  \
  --context_config single_cell/tests/jenkins/context_config.yaml \
  --submit local --loglevel DEBUG \
  --tmpdir SPLIT_WGS_BAM/temp \
  --pipelinedir SPLIT_WGS_BAM/pipeline \
  --submit local \
  --out_dir SPLIT_WGS_BAM/output \
  --config_override '{ "version": '\"$TAG\"'}'

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v $DOCKER:$DOCKER --rm \
  $1/single_cell_pipeline:$TAG \
  python single_cell/tests/jenkins/split_wgs_bam/test_split_wgs_bam.py SPLIT_WGS_BAM/output SPLIT_WGS_BAM/ref_test_data/refdata

docker run -w $PWD -v $PWD:$PWD --rm $1/single_cell_pipeline:$TAG rm -rf SPLIT_WGS_BAM
