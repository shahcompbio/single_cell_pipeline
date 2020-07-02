#!/bin/bash
set -e
set -o pipefail

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`
DOCKER=`which docker`

mkdir -p HMMCOPY/ref_test_data

docker run -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY -e AWS_DEFAULT_REGION -v $PWD:$PWD -w $PWD $1/awscli:v0.0.1 \
  aws s3 cp s3://singlecelltestsets/TESTDATA_CODEBUILD/hmmcopy HMMCOPY/ref_test_data --recursive

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v $DOCKER:$DOCKER --rm \
  $1/single_cell_pipeline:$TAG \
  single_cell hmmcopy \
  --input_yaml single_cell/tests/jenkins/hmmcopy/inputs.yaml \
  --library_id A97318A \
  --maxjobs 4 \
  --nocleanup \
  --sentinel_only  \
  --context_config single_cell/tests/jenkins/context_config.yaml \
  --submit local \
  --loglevel DEBUG \
  --config_override '{"hmmcopy": {"chromosomes": ["6", "8", "17"]}, "version": '\"$TAG\"'}' \
  --tmpdir HMMCOPY/temp \
  --pipelinedir HMMCOPY/pipeline \
  --submit local \
  --out_dir HMMCOPY/output

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v $DOCKER:$DOCKER --rm \
  $1/single_cell_pipeline:$TAG \
  python single_cell/tests/jenkins/hmmcopy/test_hmmcopy.py HMMCOPY/output A97318A  HMMCOPY/ref_test_data/refdata

docker run -w $PWD -v $PWD:$PWD --rm $1/single_cell_pipeline:$TAG rm -rf HMMCOPY
