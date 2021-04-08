#!/bin/bash
set -e
set -o pipefail

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`
DOCKER=`which docker`
NUMCORES=`nproc --all`

mkdir -p SNV_GENOTYPING/testdata

docker run -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY -e AWS_DEFAULT_REGION -v $PWD:$PWD -w $PWD $1/awscli:v0.0.1 \
  aws s3 cp s3://singlecelltestsets/TESTDATA_CODEBUILD/snv_genotyping SNV_GENOTYPING/testdata/ --recursive --quiet

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata --rm \
  $1/single_cell_pipeline_variant:$TAG \
  single_cell snv_genotyping --input_yaml single_cell/tests/jenkins/snv_genotyping/inputs.yaml \
  --maxjobs $NUMCORES --nocleanup --sentinel_only  \
  --submit local --loglevel DEBUG \
  --tmpdir SNV_GENOTYPING/temp \
  --pipelinedir SNV_GENOTYPING/pipeline --submit local --out_dir SNV_GENOTYPING/output \
  --config_override '{"variant_calling": {"chromosomes": ["6", "8", "17"]}, "version": '\"$TAG\"'}'

docker run -w $PWD -v $PWD:$PWD --rm $1/single_cell_pipeline_variant:$TAG rm -rf SNV_GENOTYPING
