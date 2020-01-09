#!/bin/bash

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`


mkdir -p VARIANT_CALLING/bams

docker run -v $PWD:$PWD -w $PWD singlecellpipeline/azurecli:v0.0.1 \
  az storage blob download-batch -s variant-calling  -d VARIANT_CALLING/bams/ --account-name $1 --account-key $2


docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v /usr/bin/docker:/usr/bin/docker --rm \
  singlecellpipeline/single_cell_pipeline:$TAG \
  single_cell variant_calling --input_yaml tests/jenkins/variant_calling/inputs.yaml \
  --maxjobs 4 --nocleanup --sentinel_only  \
  --context_config tests/jenkins/variant_calling/context_config.yaml \
  --submit local --loglevel DEBUG \
  --tmpdir VARIANT_CALLING/temp \
  --pipelinedir VARIANT_CALLING/pipeline --submit local --out_dir VARIANT_CALLING/output \
  --config_override '{"variant_calling": {"chromosomes": ["6", "8", "17"]}}'
