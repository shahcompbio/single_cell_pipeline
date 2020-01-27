#!/bin/bash

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`

#make test data directory
mkdir -p HMMCOPY/ref_test_data

#pull test data and ref data
docker run -v $PWD:$PWD -w $PWD singlecellpipeline/azurecli:v0.0.1 \
  az storage blob download-batch -s hmmcopy -d HMMCOPY/ref_test_data --account-name $1 --account-key $2

#run hmmcopy
docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v /usr/bin/docker:/usr/bin/docker --rm \
  singlecellpipeline/single_cell_pipeline:$TAG \
  single_cell hmmcopy --input_yaml tests/jenkins/hmmcopy/inputs.yaml \
  --library_id A97318A --maxjobs 4 --nocleanup --sentinel_only  \
  --context_config tests/jenkins/context_config.yaml \
  --submit local --loglevel DEBUG \
  --config_override '{"hmmcopy": {"chromosomes": ["6", "8", "17"]}}' --tmpdir HMMCOPY/temp \
  --pipelinedir HMMCOPY/pipeline --submit local --out_dir HMMCOPY/output

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v /usr/bin/docker:/usr/bin/docker --rm \
  singlecellpipeline/single_cell_pipeline:$TAG \
  python tests/jenkins/hmmcopy/test_hmmcopy.py HMMCOPY/output A97318A  HMMCOPY/ref_test_data/refdata

rm -r HMMCOPY
