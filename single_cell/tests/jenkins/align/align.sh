#!/bin/bash
set -e
set -o pipefail

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`

DOCKER=`which docker`

mkdir -p ALIGN/ref_test_data

docker run -v $PWD:$PWD -w $PWD $3/azurecli:v0.0.1 \
  az storage blob download-batch -s alignment -d ALIGN/ref_test_data --account-name $1 --account-key $2

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v $DOCKER:$DOCKER --rm \
  $3/single_cell_pipeline:$TAG \
  single_cell alignment --input_yaml single_cell/tests/jenkins/align/inputs.yaml \
  --library_id A97318A --maxjobs 4 --nocleanup --sentinel_only  \
  --context_config single_cell/tests/jenkins/context_config.yaml \
  --submit local --loglevel DEBUG \
  --tmpdir ALIGN/temp \
  --pipelinedir ALIGN/pipeline \
  --submit local \
  --out_dir ALIGN/output \
  --bams_dir ALIGN/bams \
  --config_override '{"version": '\"$TAG\"'}'


docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v $DOCKER:$DOCKER --rm \
  $3/single_cell_pipeline:$TAG \
  python single_cell/tests/jenkins/align/test_alignment.py ALIGN/output A97318A  ALIGN/ref_test_data/refdata


docker run -w $PWD -v $PWD:$PWD --rm $3/single_cell_pipeline:$TAG rm -rf ALIGN
