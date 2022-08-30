#!/bin/bash
set -e
set -o pipefail

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`
DOCKER=`which docker`
NUMCORES=`nproc --all`

mkdir -p SPLIT_WGS_BAM/ref_test_data

docker run -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY -e AWS_DEFAULT_REGION -v $PWD:$PWD -w $PWD $1/awscli:v0.0.1 \
  aws s3 cp s3://singlecelltestsets/TESTDATA_CODEBUILD/split-bam SPLIT_WGS_BAM/ref_test_data --recursive --quiet

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata --rm \
  $1/single_cell_pipeline_alignment:$TAG \
  single_cell split_wgs_bam --input_yaml single_cell/tests/codebuild/split_wgs_bam/inputs.yaml \
  --maxjobs $NUMCORES --nocleanup --sentinel_only  \
  --submit local --loglevel DEBUG \
  --tmpdir SPLIT_WGS_BAM/temp \
  --pipelinedir SPLIT_WGS_BAM/pipeline \
  --submit local \
  --output_prefix SPLIT_WGS_BAM/output/ --config_override '{"split_bam": {"chromosomes": ["6", "8", "17"]}}'

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata --rm \
  $1/single_cell_pipeline_alignment:$TAG \
  python single_cell/tests/codebuild/split_wgs_bam/test_split_wgs_bam.py SPLIT_WGS_BAM/output SPLIT_WGS_BAM/ref_test_data/refdata

docker run -w $PWD -v $PWD:$PWD --rm $1/single_cell_pipeline_alignment:$TAG rm -rf SPLIT_WGS_BAM
