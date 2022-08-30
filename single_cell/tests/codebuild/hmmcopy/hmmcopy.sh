#!/bin/bash
set -e
set -o pipefail

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`
DOCKER=`which docker`
NUMCORES=`nproc --all`

mkdir -p HMMCOPY/ref_test_data

docker run -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY -e AWS_DEFAULT_REGION -v $PWD:$PWD -w $PWD $1/awscli:v0.0.1 \
  aws s3 cp s3://singlecelltestsets/TESTDATA_CODEBUILD/hmmcopy HMMCOPY/ref_test_data --recursive --quiet

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata --rm \
  $1/single_cell_pipeline_hmmcopy:$TAG \
  single_cell hmmcopy \
  --input_yaml single_cell/tests/codebuild/hmmcopy/inputs.yaml \
  --library_id A97318A \
  --maxjobs $NUMCORES \
  --nocleanup \
  --sentinel_only  \
  --submit local \
  --loglevel DEBUG \
  --config_override '{"hmmcopy": {"chromosomes": ["6", "8", "17"]}}' \
  --tmpdir HMMCOPY/temp \
  --pipelinedir HMMCOPY/pipeline \
  --submit local \
  --output_prefix HMMCOPY/output/A97318A

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata --rm \
  $1/single_cell_pipeline_hmmcopy:$TAG \
  python single_cell/tests/codebuild/hmmcopy/test_hmmcopy.py HMMCOPY/output A97318A  HMMCOPY/ref_test_data/refdata

docker run -w $PWD -v $PWD:$PWD --rm $1/single_cell_pipeline_hmmcopy:$TAG rm -rf HMMCOPY
