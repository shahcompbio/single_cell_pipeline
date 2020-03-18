#!/bin/bash

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`


mkdir -p HMMCOPY/bams/

docker run -v $PWD:$PWD -w $PWD singlecellpipeline/azurecli:v0.0.1 \
  az storage blob download-batch -s hmmcopy -d HMMCOPY/bams/ --account-name $1 --account-key $2


docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v /usr/bin/docker:/usr/bin/docker --rm \
  singlecellpipeline/single_cell_pipeline:$TAG \
  single_cell hmmcopy --input_yaml tests/jenkins/hmmcopy/inputs.yaml \
  --library_id A97318A --maxjobs 4 --nocleanup --sentinel_only  \
  --context_config tests/jenkins/align/context_config.yaml \
  --submit local --loglevel DEBUG \
  --config_override '{"hmmcopy": {"chromosomes": ["6", "8", "17"]}}' --tmpdir HMMCOPY/temp \
  --pipelinedir HMMCOPY/pipeline --submit local --out_dir HMMCOPY/output
