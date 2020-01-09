#!/bin/bash

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`


mkdir -p MERGE_CELL_BAMS/bams/

docker run -v $PWD:$PWD -w $PWD singlecellpipeline/azurecli:v0.0.1 \
  az storage blob download-batch -s merge-bams -d MERGE_CELL_BAMS/bams/ --account-name $1 --account-key $2


docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v /usr/bin/docker:/usr/bin/docker --rm \
  singlecellpipeline/single_cell_pipeline:$TAG \
  single_cell merge_cell_bams --input_yaml tests/jenkins/merge_cell_bams/inputs.yaml \
  --maxjobs 4 --nocleanup --sentinel_only  \
  --context_config tests/jenkins/merge_cell_bams/context_config.yaml \
  --submit local --loglevel DEBUG \
  --tmpdir MERGE_CELL_BAMS/temp \
  --pipelinedir MERGE_CELL_BAMS/pipeline --submit local --out_dir MERGE_CELL_BAMS/output
