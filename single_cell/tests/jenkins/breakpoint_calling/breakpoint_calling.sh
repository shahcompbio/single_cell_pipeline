#!/bin/bash
set -e
set -o pipefail

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`
DOCKER=`which docker`
NUMCORES=`nproc --all`

mkdir -p BREAKPOINT_CALLING/ref_test_data

docker run -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY -e AWS_DEFAULT_REGION -v $PWD:$PWD -w $PWD $1/awscli:v0.0.1 \
  aws s3 cp s3://singlecelltestsets/TESTDATA_CODEBUILD/breakpoint-calling BREAKPOINT_CALLING/ref_test_data --recursive

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v $DOCKER:$DOCKER --rm \
  $1/single_cell_pipeline_breakpoint:$TAG \
  single_cell breakpoint_calling \
  --input_yaml single_cell/tests/jenkins/breakpoint_calling/inputs.yaml \
  --maxjobs $NUMCORES \
  --nocleanup \
  --sentinel_only \
  --submit local \
  --loglevel DEBUG \
  --tmpdir BREAKPOINT_CALLING/temp \
  --pipelinedir BREAKPOINT_CALLING/pipeline \
  --submit local \
  --out_dir BREAKPOINT_CALLING/output \
  --config_override '{"variant_calling": {"chromosomes": ["6", "8", "17"]}}'

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v $DOCKER:$DOCKER --rm \
  $1/single_cell_pipeline_breakpoint:$TAG \
  python single_cell/tests/jenkins/breakpoint_calling/test_breakpoint_calling.py BREAKPOINT_CALLING/output BREAKPOINT_CALLING/ref_test_data/refdata

docker run -w $PWD -v $PWD:$PWD --rm $1/single_cell_pipeline_breakpoint:$TAG rm -rf BREAKPOINT_CALLING
