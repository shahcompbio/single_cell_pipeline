#!/bin/bash
set -e
set -o pipefail

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`
TAG="${TAG}.beta"
DOCKER=`which docker`

mkdir -p VARIANT_CALLING/ref_test_data

docker run -v $PWD:$PWD -w $PWD $3/azurecli:v0.0.1 \
  az storage blob download-batch -s variant-calling  -d VARIANT_CALLING/ref_test_data/ --account-name $1 --account-key $2

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v $DOCKER:$DOCKER --rm \
  $3/single_cell_pipeline:$TAG \
  single_cell variant_calling --input_yaml single_cell/tests/jenkins/variant_calling/inputs.yaml \
  --maxjobs 4 --nocleanup --sentinel_only  \
  --context_config single_cell/tests/jenkins/context_config.yaml \
  --submit local --loglevel DEBUG \
  --tmpdir VARIANT_CALLING/temp \
  --pipelinedir VARIANT_CALLING/pipeline --submit local --out_dir VARIANT_CALLING/output \
  --config_override '{"variant_calling": {"chromosomes": ["6", "8", "17"]}, "version": '\"$TAG\"'}'

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v $DOCKER:$DOCKER --rm \
  $3/single_cell_pipeline:$TAG \
  python single_cell/tests/jenkins/variant_calling/test_variant_calling.py VARIANT_CALLING/output VARIANT_CALLING/ref_test_data/refdata

docker run -w $PWD -v $PWD:$PWD --rm $3/single_cell_pipeline:$TAG rm -rf VARIANT_CALLING
