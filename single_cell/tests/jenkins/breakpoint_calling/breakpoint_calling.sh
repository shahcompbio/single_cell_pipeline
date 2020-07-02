#!/bin/bash
set -e
set -o pipefail

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`
DOCKER=`which docker`

mkdir -p BREAKPOINT_CALLING/ref_test_data

docker run -v $PWD:$PWD -w $PWD $3/azurecli:v0.0.1 \
  az storage blob download-batch -s breakpoint-calling  -d BREAKPOINT_CALLING/ref_test_data --account-name $1 --account-key $2

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v $DOCKER:$DOCKER --rm \
  $3/single_cell_pipeline:$TAG \
  single_cell breakpoint_calling \
  --input_yaml single_cell/tests/jenkins/breakpoint_calling/inputs.yaml \
  --maxjobs 4 \
  --nocleanup \
  --sentinel_only \
  --context_config single_cell/tests/jenkins/context_config.yaml \
  --submit local \
  --loglevel DEBUG \
  --tmpdir BREAKPOINT_CALLING/temp \
  --pipelinedir BREAKPOINT_CALLING/pipeline \
  --submit local \
  --out_dir BREAKPOINT_CALLING/output \
  --config_override '{"variant_calling": {"chromosomes": ["6", "8", "17"]}, "version": '\"$TAG\"'}'

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v $DOCKER:$DOCKER --rm \
  $3/single_cell_pipeline:$TAG \
  python single_cell/tests/jenkins/breakpoint_calling/test_breakpoint_calling.py BREAKPOINT_CALLING/output BREAKPOINT_CALLING/ref_test_data/refdata

docker run -w $PWD -v $PWD:$PWD --rm $3/single_cell_pipeline:$TAG rm -rf BREAKPOINT_CALLING
